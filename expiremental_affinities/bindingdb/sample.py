"""
sample.py
Query BindingDB for protein-ligand complexes similar to a query ligand,
randomly sample N of them, fetch protein sequences from UniProt.
"""


from rdkit import Chem, DataStructs
from concurrent.futures import ProcessPoolExecutor
import os
import sys
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from rdkit.Chem import rdFingerprintGenerator
from chemspipy import ChemSpider

mfp = rdFingerprintGenerator.GetMorganGenerator()
CS = ChemSpider("TzfV80skfR1yCAe7oY4y06Q2s5Kzlff4a1NTNVq2")
BINDINGDB_URL = "https://bindingdb.org/rest/getTargetByCompound"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{uid}.fasta"
INTRA_LIGAND_WORKERS = 16
INTER_LIGAND_WORKERS = 8

def smiles_from_chemspider(chemspi_id: str) -> str:   
    smiles = CS.get_details(chemspi_id)['smiles']
    return smiles


def smiles_from_csv(csv_path: Path | str) -> dict[str, str]:

    df = pd.read_csv(csv_path)

    cols = {c.lower(): c for c in df.columns}
    chem_col = cols.get("chemspider id")
    smiles_col = cols.get("smiles")
    name_col = cols.get("common name")

    if smiles_col not in df.columns:
        df[smiles_col] = ""


    all_smiles: dict[str, str] = {}
    changed = False

    for idx, row in df.iterrows():
        common_name = str(row.get(name_col)) 
        chem_id = str(row.get(chem_col, "")).strip()
        smiles = str(row.get(smiles_col, "")).strip()

        if not smiles:
            smiles = smiles_from_chemspider(chem_id) or ""
            if smiles:
                df.at[idx, smiles_col] = smiles
                changed = True

        all_smiles[f"{common_name}_{chem_id}"] = smiles

    if changed:
        df.to_csv(csv_path, index=False)

    return all_smiles


_WORKER_Q_FP = None


def _init_similarity_worker(query_smiles: str):
    global _WORKER_Q_FP
    _WORKER_Q_FP = _morgan_fp(query_smiles)
    if _WORKER_Q_FP is None:
        raise ValueError(f"could not parse query SMILES: {query_smiles!r}")


def _score_smiles_tanimoto(smi: str) -> tuple[str, float]:
    fp = _morgan_fp(smi)
    if fp is None:
        return smi, 0.0
    return smi, DataStructs.TanimotoSimilarity(_WORKER_Q_FP, fp)  # type: ignore


def _morgan_fp(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return mfp.GetFingerprint(mol)


def sample_complexes(
    index_path: Path | str,
    query_smiles: str,
    n: int,
    similarity_cutoff: float,
    seed: int = 42,
    require_affinity: bool = True,
    organism_contains: str | None = None,
    fingerprint_workers: int = INTER_LIGAND_WORKERS,
    show_progress: bool = True,
) -> pd.DataFrame:
    df = pd.read_parquet(index_path)
    print(f"[info] index: {len(df):,} measurements", file=sys.stderr)

    if require_affinity:
        aff = ["ki_nM", "ic50_nM", "kd_nM", "ec50_nM"]
        df = df[df[aff].notna().any(axis=1)]
    if organism_contains:
        df = df[df["organism"].fillna("").str.contains(organism_contains, case=False)]
    print(f"[info] after filters: {len(df):,}", file=sys.stderr)
    if df.empty:
        return df

    q_fp = _morgan_fp(query_smiles)
    if q_fp is None:
        raise ValueError(f"could not parse query SMILES: {query_smiles!r}")

    unique_smiles = df["smiles"].drop_duplicates().tolist()
    print(
        f"[info] fingerprinting {len(unique_smiles):,} unique ligands", file=sys.stderr
    )

    workers = max(1,min(fingerprint_workers,len(unique_smiles)))

    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=_init_similarity_worker,
        initargs=(query_smiles,),
    ) as ex:
        sims = dict(
            tqdm(
                ex.map(_score_smiles_tanimoto, unique_smiles, chunksize=500),
                total=len(unique_smiles),
                desc="fingerprints",
                unit="mol",
                disable=not show_progress,
            )
        )

    df = df.assign(tanimoto=df["smiles"].map(sims))
    hits = df[df["tanimoto"] >= similarity_cutoff].copy()
    print(
        f"[info] {len(hits):,} measurements >= cutoff {similarity_cutoff}",
        file=sys.stderr,
    )
    if hits.empty:
        return hits

    hits["target_key"] = hits["uniprots"].map(
        lambda xs: "|".join([u or "" for u in xs])  # type: ignore
    )
    hits = (
        hits.sort_values("tanimoto", ascending=False)
        .drop_duplicates(subset=["target_key", "monomer_id"], keep="first")
        .drop(columns=["target_key"])
    )
    print(f"[info] {len(hits):,} unique (target, ligand) complexes", file=sys.stderr)

    if len(hits) > n:
        hits = hits.sample(n=n, random_state=seed)
    return hits.reset_index(drop=True)


def _run_ligand_query(args: tuple[str, str]) -> tuple[str, int, pd.DataFrame]:
    name, smi = args

    hits = sample_complexes(
        "bindingdb_index.parquet",
        query_smiles=smi,
        n=10,
        similarity_cutoff=1.0,
        require_affinity=True,
        fingerprint_workers=INTRA_LIGAND_WORKERS,
        show_progress=False,
    )

    return name, len(hits), hits


if __name__ == "__main__":
    ligands = smiles_from_csv("ligands.csv")
    print(f"got {len(ligands)} ligands from CSV")

    total_cores = os.cpu_count() or 8
    inter_workers = min(INTER_LIGAND_WORKERS, len(ligands))
    INTRA_LIGAND_WORKERS = max(1, total_cores // max(1, inter_workers))

    print(
        f"using {inter_workers} ligand workers x "
        f"{INTRA_LIGAND_WORKERS} fingerprint workers each"
    )

    with ProcessPoolExecutor(max_workers=inter_workers) as ex:
        results = list(
            tqdm(
                ex.map(_run_ligand_query, ligands.items()),
                total=len(ligands),
                desc="ligands",
                unit="ligand",
            )
        )

    for name, n_hits, hits in results:
        print(f"found {n_hits} hits for {name}")
