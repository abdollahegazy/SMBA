import pandas as pd
from pathlib import Path
import subprocess

from utils import DATA_DIR

def extract_best_ligand(out_poses:Path):
    content = out_poses.read_text()
    models = {}

    current_model = None
    current_affinity = None
    current_lines = []

    for line in content.splitlines():
        if line.startswith("MODEL"):
            current_model = int(line.split()[1])
            current_lines = [line]
        elif line.startswith("REMARK VINA RESULT:"):
            current_affinity = float(line.split()[3])
            current_lines.append(line)
        elif line.startswith("ENDMDL"):
            current_lines.append(line)
            models[current_model] = {
                "ligand_file": "\n".join(current_lines),
                "affinity": current_affinity,
            }
        else:
            current_lines.append(line)

    return models

def create_complex(protein_path:Path,ligand_path:Path,outpath:Path):
    subprocess.run(["vmd2","-dispdev","text","-e","smba_post/combine.tcl","-args",protein_path,ligand_path,outpath],check=True)


def main():
    data = pd.read_csv(DATA_DIR / 'docking_results.csv')

    for _, row in data.iterrows():
        protein = row['protein']
        ligand = row['ligand']
        affinity = row["best_affinity"]
        conformation = row['conformation']
        pocket = row['pocket']

        protein_dir = (DATA_DIR/protein).resolve()
        structures_dir = Path(f'{protein_dir}/smba_structures')
        structures_dir.mkdir(exist_ok=True,parents=True)

        outpath = structures_dir /  f"{ligand}_c{conformation}_p{pocket}.pdb" 
        
        # if exists and more than 10 lines. kinda arbitrary
        if outpath.exists() and len(outpath.read_text().splitlines()) > 10:
            continue
        
        smba_ligands_path = protein_dir / "docking" / str(ligand) / f"out_poses_c{conformation}_p{pocket}.pdbqt"
        ligand_models = extract_best_ligand(smba_ligands_path)
        best_ligand,best_affiinity = ligand_models[1].values()
        assert best_affiinity == affinity, "Something went wrong extracting the best ligand"
    
        ligand_filename = f"{ligand}_c{conformation}_p{pocket}.pdbqt"
        ligand_filepath = structures_dir / ligand_filename
        ligand_filepath.write_text(best_ligand)

        corresponding_protein = protein_dir / "conformation" / f"protein_conf{conformation}.pdb"

        create_complex(corresponding_protein,ligand_filepath,outpath)
        ligand_filepath.unlink() 




if __name__ == "__main__":
    main()