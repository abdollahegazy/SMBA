"""Prepare per-conformation pdbqts for the AF3-pocket experiment.

For each SMBA conformation i, RMSD-aligns the AF3 protein+ligand structure onto
protein_conf{i}.pdb, takes the aligned AF3 ligand COM as the conf-specific pocket
center, and writes a Meeko-prepped pdbqt pre-centered at that COM:

    ligands_prepared_af3_pocket/{lid}/ligand_c{i}_p0.pdbqt

Also caches the AF3 protein+ligand structure as PDB at
    inference_af3_lig/{lid}/{lid}.pdb
so this script (and downstream consumers) can read it with MDAnalysis without
re-converting.
"""

from argparse import ArgumentParser, Namespace
from pathlib import Path

import gemmi
from MDAnalysis import Universe

from smba.prepare_ligand import smiles_to_pdb
from varidock.stages.center_ligand_to_pocket import CenterLigand, CenterLigandConfig
from varidock.stages.meeko_ligand_prep import MeekoLigandPrep, MeekoLigandPrepConfig
from varidock.types import Ligand, LigandPrepInput, PDB, PocketCenter

from utils import DATA_DIR, load_pairs, ligand_center, rmsd_align

NUM_CONFS = 11


def load_args() -> Namespace:
    parser = ArgumentParser(
        description="Generate per-conf pdbqts (pre-centered at AF3 ligand COM) for AF3-pocket redocking."
    )
    return parser.parse_args()


def af3_protein_ligand_pdb(protein_id: str, ligand_id: str) -> Path:
    """Convert AF3 protein+ligand CIF to PDB (cached). Preserves chain L/P."""
    inf_dir = DATA_DIR / protein_id / "inference_af3_lig" / ligand_id
    pdb_path = inf_dir / f"{ligand_id}.pdb"
    if pdb_path.exists() and pdb_path.stat().st_size > 0:
        return pdb_path

    cif_path = inf_dir / "af_output" / ligand_id / f"{ligand_id}_model.cif"
    if not cif_path.exists():
        raise FileNotFoundError(f"AF3 CIF not found: {cif_path}")

    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    st.write_pdb(str(pdb_path))
    return pdb_path


def af3_ligand_com_in_conf_frame(
    af3_pdb: Path, conf_pdb: Path
) -> tuple[float, float, float]:
    """RMSD-align AF3 protein onto conf, return AF3 ligand COM in conf frame."""
    conf_u = Universe(str(conf_pdb.resolve()))
    af3_u = Universe(str(af3_pdb.resolve()))
    rmsd_align(
        moving_universe=af3_u,
        reference_universe=conf_u,
        moving_sel="chainID P and backbone and name CA",
        reference_sel="backbone and name CA",
    )
    return ligand_center(af3_u, "chainID L")


def run_one(protein_id: str, ligand_id: str, smiles: str) -> None:
    ligand_dir = DATA_DIR / protein_id / "ligands_prepared_af3_pocket" / ligand_id
    expected_pdbqts = [
        ligand_dir / f"ligand_c{i}_p0.pdbqt" for i in range(NUM_CONFS)
    ]
    if all(p.exists() and p.stat().st_size > 0 for p in expected_pdbqts):
        return

    af3_pdb = af3_protein_ligand_pdb(protein_id, ligand_id)

    ligand_dir.mkdir(parents=True, exist_ok=True)
    raw_pdb = ligand_dir / f"{ligand_id}.pdb"
    if not raw_pdb.exists():
        smiles_to_pdb(smiles, raw_pdb)

    for conf_idx in range(NUM_CONFS):
        pdbqt_path = ligand_dir / f"ligand_c{conf_idx}_p0.pdbqt"
        if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
            continue

        conf_pdb = (
            DATA_DIR / protein_id / "conformation" / f"protein_conf{conf_idx}.pdb"
        )
        if not conf_pdb.exists():
            continue

        com = af3_ligand_com_in_conf_frame(af3_pdb, conf_pdb)

        inp = LigandPrepInput(
            ligand=Ligand(name=ligand_id, pdb=PDB(path=raw_pdb.resolve())),
            pocket_center=PocketCenter(
                x=float(com[0]), y=float(com[1]), z=float(com[2])
            ),
            conf_index=conf_idx,
            pose_index=0,
        )
        centered = CenterLigand(
            CenterLigandConfig(output_dir=ligand_dir.resolve())
        ).run(inp)
        MeekoLigandPrep(
            MeekoLigandPrepConfig(output_dir=ligand_dir.resolve())
        ).run(centered)

        # MeekoLigandPrep leaves a *_protonated.pdb intermediate, and CenterLigand
        # leaves the centered .pdb -- drop both, keep only the .pdbqt.
        (ligand_dir / f"ligand_c{conf_idx}_p0.pdb").unlink(missing_ok=True)
        (ligand_dir / f"ligand_c{conf_idx}_p0_protonated.pdb").unlink(missing_ok=True)


def main():
    load_args()
    pairs = load_pairs()

    for protein_id, ligand_id, smiles in pairs:
        run_one(protein_id, ligand_id, smiles)


if __name__ == "__main__":
    main()
