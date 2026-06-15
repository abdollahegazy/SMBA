from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import AllChem
from utils import load_pairs, DATA_DIR
import traceback

from varidock.stages.meeko_ligand_prep import MeekoLigandPrep, MeekoLigandPrepConfig
from varidock.stages.center_ligand_to_pocket import CenterLigand, CenterLigandConfig

from varidock.types import LigandPrepInput, PocketCenter, Ligand, PDB
MAX_WORKERS = 64
NUM_CONFS = 11



# def smiles_to_pdb(smiles: str, output_path: Path) -> None:
#     mol = Chem.MolFromSmiles(smiles)
#     if mol is None:
#         raise ValueError(f"Invalid SMILES: {smiles}")
#     mol = Chem.AddHs(mol)
#     AllChem.EmbedMolecule(mol, randomSeed=42) #type: ignore
#     AllChem.MMFFOptimizeMolecule(mol) #type: ignore
#     Chem.MolToPDBFile(mol, str(output_path))

# atropine is rlly weird and would fail for this (also did smth similar for boltz)
# so we are doing this hack
def smiles_to_pdb(smiles: str, output_path: Path) -> None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)

    # Snapshot stereo for re-imposition if we fall back
    original_atom_chirality = [a.GetChiralTag() for a in mol.GetAtoms()]
    original_bond_stereo = [b.GetStereo() for b in mol.GetBonds()]

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    conf_id = AllChem.EmbedMolecule(mol, params)

    if conf_id == -1:
        params.useRandomCoords = True
        conf_id = AllChem.EmbedMolecule(mol, params)

    if conf_id == -1:
        params.enforceChirality = False
        params.ignoreSmoothingFailures = True
        conf_id = AllChem.EmbedMolecule(mol, params)
        if conf_id != -1:
            # re-impose stereo on the relaxed embed
            for atom, tag in zip(mol.GetAtoms(), original_atom_chirality):
                atom.SetChiralTag(tag)
            for bond, stereo in zip(mol.GetBonds(), original_bond_stereo):
                bond.SetStereo(stereo)

    if conf_id == -1:
        raise ValueError(f"Failed to embed conformer for SMILES: {smiles}")

    # Try MMFF; fall back to UFF since MMFF doesn't cover all atom types
    try:
        if AllChem.MMFFOptimizeMolecule(mol, confId=conf_id, maxIters=2000) == -1:
            AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=2000)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=2000)

    Chem.MolToPDBFile(mol, str(output_path))



def process_protein(protein_id: str, ligand_id: str, smiles: str) -> None:
    protein_dir = DATA_DIR / protein_id
    pocket_dir = protein_dir / "pockets"
    ligand_dir = protein_dir / "ligands_prepared" / ligand_id
    ligand_dir.mkdir(parents=True, exist_ok=True)

    # Generate ligand PDB from SMILES
    ligand_pdb = ligand_dir / f"{ligand_id}.pdb"
    if not ligand_pdb.exists():
        smiles_to_pdb(smiles, ligand_pdb)

    for i in range(NUM_CONFS):
        centers_file = pocket_dir / f"protein_conf{i}" / "centers.txt"
        if not centers_file.exists() or centers_file.stat().st_size == 0:
            continue

        with open(centers_file) as f:
            for pocket, line in enumerate(f):
                x, y, z = map(float, line.strip().split())

                pdbqt = ligand_dir / f"ligand_c{i}_p{pocket}.pdbqt"
                if pdbqt.exists() and pdbqt.stat().st_size > 0:
                    continue

                inp = LigandPrepInput(
                    ligand=Ligand(name=ligand_id, pdb=PDB(path=ligand_pdb.resolve())),
                    pocket_center=PocketCenter(x=x, y=y, z=z),
                    conf_index=i,
                    pose_index=pocket,
                )

                center_config = CenterLigandConfig(output_dir=ligand_dir.resolve())
                centered = CenterLigand(center_config).run(inp)

                meeko_config = MeekoLigandPrepConfig(output_dir=ligand_dir.resolve())
                MeekoLigandPrep(meeko_config).run(centered)

                # cleanup intermediates
                centered_pdb = ligand_dir / f"ligand_c{i}_p{pocket}.pdb"
                centered_pdb.unlink(missing_ok=True)
                protonated = ligand_dir / f"ligand_c{i}_p{pocket}_protonated.pdb"
                protonated.unlink(missing_ok=True)


def main():
    pairs = load_pairs()

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = [
            pool.submit(process_protein, pid, lid, smi) for pid, lid, smi in pairs
        ]

        for fut in tqdm(as_completed(futures), total=len(futures), desc="Ligand prep"):
            try:
                fut.result()
            except Exception as e:
                print(f"ERROR: {e}")
                traceback.print_exc()


if __name__ == "__main__":
    main()
