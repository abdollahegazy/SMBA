
from __future__ import annotations

# from concurrent.futures import ProcessPoolExecutor
import csv
# import os
import sys
from pathlib import Path
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm
import re
# from rdkit import Chem


CORE_COLS = {
    "BindingDB Reactant_set_id": "measurement_id",
    "BindingDB MonomerID": "monomer_id",
    "Ligand SMILES": "smiles",
    "Target Name": "target_name",
    "Target Source Organism According to Curator or DataSource": "organism",
    "Ki (nM)": "ki_nM",
    "IC50 (nM)": "ic50_nM",
    "Kd (nM)": "kd_nM",
    "EC50 (nM)": "ec50_nM",
    "pH": "pH",
    "Temp (C)": "temp_C",
    "PDB ID(s) for Ligand-Target Complex": "pdb_complex",
    "Number of Protein Chains in Target (>1 implies a multichain complex)": "num_chains",
    "Ligand InChI Key": "inchikey",
}

# chain-N column name templates
CHAIN_SEQ_RE = re.compile(r"^BindingDB Target Chain Sequence (\d+)$")
CHAIN_UNIPROT_RE = re.compile(
    r"^UniProt \(SwissProt\) Primary ID of Target Chain (\d+)$"
)
CHAIN_PDB_RE = re.compile(r"^PDB ID\(s\) of Target Chain (\d+)$")

PARQUET_SCHEMA = pa.schema(
    [
        ("measurement_id", pa.string()),
        ("monomer_id", pa.string()),
        ("smiles", pa.string()),
        ("target_name", pa.string()),
        ("organism", pa.string()),
        ("ki_nM", pa.string()),
        ("ic50_nM", pa.string()),
        ("kd_nM", pa.string()),
        ("ec50_nM", pa.string()),
        ("pH", pa.string()),
        ("temp_C", pa.string()),
        ("pdb_complex", pa.string()),
        ("num_chains", pa.string()),
        ("sequences", pa.list_(pa.string())),
        ("uniprots", pa.list_(pa.string())),
        ("pdb_chains", pa.list_(pa.string())),
        ("inchikey", pa.string()),
    ]
)



# ---------------- index build ----------------


def _discover_chain_cols(header: list[str]) -> dict[int, dict[str, str]]:
    """Return {chain_number: {'seq': col, 'uniprot': col, 'pdb': col}, ...}."""
    by_n: dict[int, dict[str, str]] = {}
    for c in header:
        for kind, regex in [
            ("seq", CHAIN_SEQ_RE),
            ("uniprot", CHAIN_UNIPROT_RE),
            ("pdb", CHAIN_PDB_RE),
        ]:
            m = regex.match(c)
            if m:
                n = int(m.group(1))
                by_n.setdefault(n, {})[kind] = c
    return by_n


def build_index(tsv_path: Path | str, out_path: Path | str, batch_size: int = 25_000) -> None:
    """ Stream the BindingDB TSV and write a Parquet to be faster (one row per measurement)."""

    if isinstance(tsv_path, str):
        tsv_path = Path(tsv_path)
    if isinstance(out_path, str):
        out_path = Path(out_path)

    with tsv_path.open(newline="", encoding="utf-8", errors="replace") as f:
        header = next(csv.reader(f, delimiter="\t"))

    chain_cols = _discover_chain_cols(header)
    
    writer = pq.ParquetWriter(str(out_path), PARQUET_SCHEMA, compression="zstd")
    batch: list[dict] = []
    n_total = n_kept = 0

    def flush(rows: list[dict]) -> None:
        if not rows:
            return
        df = pd.DataFrame(rows, columns=PARQUET_SCHEMA.names)
        writer.write_table(
            pa.Table.from_pandas(df, schema=PARQUET_SCHEMA, preserve_index=False)
        )

    with tsv_path.open(newline="", encoding="utf-8", errors="replace") as f:
        # csv module handles quoting/embedded tabs better than naive split
        reader = csv.DictReader(f, delimiter="\t")
        for raw in tqdm(reader, desc="parsing TSV", unit="row", total=3_176_529):
            # skip row if num_chains is more than 1
            
            if raw.get("Number of Protein Chains in Target (>1 implies a multichain complex)", "").strip() != "1":
                continue

            n_total += 1
            smi = (raw.get("Ligand SMILES") or "").strip()
            if not smi:
                continue

            # Collect chains in order, skipping empties
            # ony get the first chain, since we skip multichain complexes
            seqs, uids, pdbs = [], [], []
            for n in sorted(chain_cols):
                cols = chain_cols[n]
                seq = (raw.get(cols.get("seq", "")) or "").strip()
                if not seq:
                    continue
                seqs.append(seq)
                uids.append((raw.get(cols.get("uniprot", "")) or "").strip() or None)
                pdbs.append((raw.get(cols.get("pdb", "")) or "").strip() or None)
                break

            if not(len(seqs) == 1 and len(uids) == 1 and len(pdbs) == 1):
                continue  # dont care ab it ngl 

            row = {out: (raw.get(src) or None) for src, out in CORE_COLS.items()}
            # strip empty strings -> None
            row = {
                k: (v.strip() if isinstance(v, str) else v) or None
                for k, v in row.items()
            }
            row["smiles"] = smi

            row["sequences"] = seqs
            row["uniprots"] = uids
            row["pdb_chains"] = pdbs
            batch.append(row)
            n_kept += 1

            if len(batch) >= batch_size:
                flush(batch)
                batch = []

    flush(batch)
    writer.close()
    print(f"[done] parsed {n_total} rows, kept {n_kept} -> {out_path}", file=sys.stderr)


build_index("BindingDB_All.tsv", "bindingdb.parquet", batch_size=50_000)

# def canonical_smiles(smi: str) -> str | None:
#     mol = Chem.MolFromSmiles(str(smi))
#     if mol is None:
#         return None
#     return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


# def canonical_pair(smi: str) -> tuple[str, str | None]:
#     return smi, canonical_smiles(smi)


# if __name__ == "__main__":
#     in_path = Path("bindingdb.parquet")
#     out_path = Path("bindingdb_with_canonical.parquet")

#     df = pd.read_parquet(in_path)
#     print(f"loaded {len(df):,} rows")

#     unique_smiles = df["smiles"].dropna().drop_duplicates().tolist()
#     print(f"canonicalizing {len(unique_smiles):,} unique SMILES")

#     workers = max(1, (os.cpu_count() or 2) - 1)

#     with ProcessPoolExecutor(max_workers=workers) as ex:
#         canon_map = dict(
#             tqdm(
#                 ex.map(canonical_pair, unique_smiles, chunksize=1000),
#                 total=len(unique_smiles),
#                 desc="canonicalizing",
#                 unit="mol",
#             )
#         )

#     df["canonical_smiles"] = df["smiles"].map(canon_map)

#     df.to_parquet(out_path, index=False, compression="zstd")
#     print(f"wrote {out_path}")