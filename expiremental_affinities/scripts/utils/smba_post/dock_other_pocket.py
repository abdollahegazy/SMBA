from pathlib import Path
from MDAnalysis import Universe
import numpy as np
from utils import DATA_DIR, rmsd_align,ligand_center
from varidock.stages.vina_dock import VinaDocking, VinaDockingConfig
from varidock.types import DockingInput, PocketCenter, PDBQT
from typing import Tuple
from argparse import ArgumentParser,Namespace

def parse_args() -> Namespace:
    parser = ArgumentParser(description="Dock ligand into SMBA pocket centered at Boltz pose")
    parser.add_argument("--protein_id", type=str, required=True, help="Protein ID")
    parser.add_argument("--ligand_id", type=str, required=True, help="Ligand ID")
    parser.add_argument("--conformation_idx", type=int, required=True, help="Conformation index")
    parser.add_argument("--pocket_idx", type=int, required=True, help="Pocket index")
    parser.add_argument("--reference", type=str, choices=["boltz"], default="boltz", help="Reference structure to align to")
    return parser.parse_args()


def get_target_ligand_center(
        protein_id: str, 
        ligand_id: str, 
        conformation_idx: int,
        pocket_idx: int,
        reference: str) -> tuple[float, float, float]:
    if reference == "boltz":
        other_path = Path(f"/serviceberry/tank/abdolla/SMBA/expiremental_affinities/data/{protein_id}/boltz/{ligand_id}/boltz_output/boltz_results_{protein_id}_{ligand_id}/predictions/{protein_id}_{ligand_id}/{protein_id}_{ligand_id}_model_0.pdb")
    else:
        other_path = Path()
    smba_path = Path(
        f"/serviceberry/tank/abdolla/SMBA/expiremental_affinities/data/{protein_id}/smba_structures/{ligand_id}_c{conformation_idx}_p{pocket_idx}.pdb"
    )
    assert other_path.exists(), f"Reference structure not found at {other_path}"
    assert smba_path.exists(), f"SMBA structure not found at {smba_path}. Maybe non-existent conformation/pocket."

    #this stays stationary
    smba_universe = Universe(str(smba_path.resolve()))
    #this is aligned onto SMBA 
    other_universe = Universe(str(other_path.resolve()))

    rmsd_align(
        moving_universe=other_universe,
        reference_universe=smba_universe,
        moving_sel="backbone and name CA",
        reference_sel="backbone and name CA"
    )

    # other ligand center
    return ligand_center(other_universe, "chainID L")

def generate_vina_ligand_at_other_center(
        protein_id: str, 
        ligand_id: str, 
        conformation_idx: int,
        pocket_idx: int,
        reference: str) -> Tuple[Path, PocketCenter]:
    target_center = np.array(get_target_ligand_center(protein_id, ligand_id, conformation_idx, pocket_idx, reference),dtype=float)

    protein_dir = DATA_DIR / protein_id
    ligand_prepared_dir = protein_dir / "ligands_prepared" / ligand_id
    ligand_pdbqt = ligand_prepared_dir / f"ligand_c{conformation_idx}_p{pocket_idx}.pdbqt"

    if not ligand_pdbqt.exists():
        raise ValueError(
            f"Missing ligand PDBQT for {ligand_id} in {protein_id} for c={conformation_idx} p={pocket_idx}"
        )
    
    # read ligand pdbqt, shift by difference in centers, write new pdbqt
    lines = ligand_pdbqt.read_text().splitlines()
    atom_idx = []
    coords = []
    for i, line in enumerate(lines):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x, y, z = map(float, [line[30:38], line[38:46], line[46:54]])
            atom_idx.append(i)
            coords.append([x, y, z])
    coords = np.array(coords, dtype=float)
    current_center = coords.mean(axis=0)  
    shift = target_center - current_center

    # apply translation
    for j, i in enumerate(atom_idx):
        x, y, z = coords[j] + shift
        # print(f"Original coords for atom {j}: {coords[j]}, shifted coords: {[x,y,z]}")
        if len(lines[i]) < 54:
            lines[i] = lines[i].rstrip("\n")
            lines[i] = lines[i] .ljust(54)
        # dont touch anything except XYZ cols
        lines[i] = f"{lines[i][:30]}{x:8.3f}{y:8.3f}{z:8.3f}{lines[i][54:]}"

    out_dir = protein_dir / "ligands_prepared" / f"{ligand_id}_{reference}_centered"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / f"ligand_c{conformation_idx}_p{pocket_idx}.pdbqt"
    out_path.write_text("\n".join(lines) + "\n")

    # print("target_center:", target_center)
    # print("old_center:   ", current_center)
    # print("new_center:   ", current_center + shift)
    # print("wrote:", out_path)
    return out_path, PocketCenter(x=target_center[0], y=target_center[1], z=target_center[2])


def docking(
    protein_id: str,
    conformation_idx: int,
    pocket_idx: int,
    ligand_path: Path,
    center: PocketCenter,
    out_dir: Path
):
    out_dir.mkdir(exist_ok=True, parents=True)
    log_file = out_dir / f"out_c{conformation_idx}_p{pocket_idx}.log"
    if log_file.exists() and log_file.stat().st_size > 0:
        print(f"Skipping docking for {protein_id} c{conformation_idx} p{pocket_idx} - log already exists")
        return  
    receptor = PDBQT(
        Path(f"../data/{protein_id}/protein_receptors/protein_conf{conformation_idx}.pdbqt")
    )
    ligand = PDBQT(path=ligand_path)

    config = VinaDockingConfig(output_dir=out_dir.resolve())
    VinaDocking(config).run(
        DockingInput(
            receptor=receptor,
            ligand=ligand,
            pocket_center=center,
            conf_index=conformation_idx,
            pose_index=pocket_idx,
        )
    )

def main():
    args = parse_args()
    protein_id = args.protein_id
    ligand_id = args.ligand_id
    conformation_idx = args.conformation_idx
    pocket_idx = args.pocket_idx
    reference = args.reference

    ligand_path, pocket_center = generate_vina_ligand_at_other_center(
        protein_id, ligand_id, conformation_idx, pocket_idx, reference)
    
    docking_dir = ligand_path.parents[2] / "docking" / ligand_path.parent.name

    docking(
        protein_id=protein_id,
        conformation_idx=conformation_idx,
        pocket_idx=pocket_idx,
        ligand_path=ligand_path,
        center=pocket_center,
        out_dir=docking_dir
    )


if __name__ == "__main__":
    main()
