from pathlib import Path
import json
from pprint import pprint 
import time
import os
from tqdm import tqdm
from utils import load_ligand_data,AF3_INPUT_DIR

OUT_DIR = Path("../input")



def create_dirs(ligand_data:dict[str:str]) -> None:
    """
    Iterates all the protein data from AF3 docking input and basically mirrors structure to create Boltz-2 input directories here.

    Args:
        ligand data (dict): mapping of ligand ids and their smiles.
    """
    all_dirs_to_create = []
    for species in AF3_INPUT_DIR.iterdir():
        for protein in species.iterdir():
            if protein.name in ['scripts','batches']:continue #remenant logging file we don't care about 
            
            for ligand_name in ligand_data.keys():
                ligand_out = OUT_DIR / species.name / protein.name / ligand_name
                all_dirs_to_create.append(ligand_out)
        #         break
        #     break
        # break

    print(f"Creating {len(all_dirs_to_create)} directories...");s = time.time()
    for dir_path in all_dirs_to_create:
        dir_path.mkdir(parents=True, exist_ok=True)
    print(f"Creating {len(all_dirs_to_create)} dirs took {time.time()-s:.2f} seconds.")


def process_msas(protein_dir:Path, unpairedmsa:str) -> None:
    """
    Writes the unpaired MSA to its respective protein in the boltz-2 input files.

    Args:
        protein_dir: Path of where the protein data is being saved to in Boltz-2 input files
    """
    with open(str(protein_dir/"unpaired_msa.a3m"),'w') as f:
        f.write(unpairedmsa)


def process_protein(protein_fasta:str,
                    ligand_out:Path,
                    ligand_smiles:str) -> None:
    """
    Writes the Boltz-2 yaml file to the respective ligand dir in the boltz-2 input files.

    Args:
        protein_fasta: string of the protein fasta to include
        ligand_out: Path of where the ligand output file is
        ligand_smiles: string of the ligand smiles to include
    """

    yaml = f"""
version: 1 
sequences:
    - protein: 
        id: A
        sequence: {protein_fasta}
        msa: ../unpaired_msa.a3m
    - ligand:
        id: B
        smiles: '{ligand_smiles}'
properties:
    - affinity:
        binder: B
"""
    
    with open(str(ligand_out / "run.yaml"),'w') as f:
        f.write(yaml.lstrip())


def process_proteins(ligand_data:dict[str:str]):
    """
    Iterates all the protein directories to create their respective yaml files. It yanks the majority of the info by loading the first AF3 json equivalent of that protein.

    Args:
        ligand_data: mapping of ligand_id : ligand_smiles
    """

    for species in AF3_INPUT_DIR.iterdir():
        protein_files = list(species.iterdir())
        for protein in tqdm(protein_files,desc=f'Creating YAML"s for  {species.name}"s'):
            
            if protein.name in ['scripts','batches']:continue #remenant logging file we don't care about 


            # readad protein data once per protein
            first_ligand_file = next(protein.glob("*.json"), None)

            if first_ligand_file is None and (protein.name!= 'scripts' or protein.name!='batches'):
                raise FileNotFoundError(f"ERROR! No JSON found", protein)

            with open(first_ligand_file) as f:
                data = json.load(f)

            protein_info = data['sequences'][0]['protein']
            
            fasta = protein_info['sequence']
            unpairedmsa = protein_info['unpairedMsa']
            
            protein_out =(OUT_DIR / species.name / protein.name) 
            process_msas(protein_out,unpairedmsa)

            for ligand_name, ligand_smiles in ligand_data.items():
                ligand_out = OUT_DIR / species.name / protein.name / ligand_name

                process_protein(fasta,ligand_out,ligand_smiles)


if __name__ == "__main__":
    ligand_mapping = load_ligand_data()
    create_dirs(ligand_mapping)
    process_proteins(ligand_mapping)