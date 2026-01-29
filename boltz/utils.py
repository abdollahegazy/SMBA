from pathlib import Path
import json

AF3_INPUT_DIR = Path("../../dockingAF3/input/")

def load_ligand_data() -> dict[str:str]:
    """
    Reads the first protein in AF3 path to scan all its dirs for all the ligands we use for testing.
    Collects the ligand ID (the json filename) and its smiles from the json

    Args:
        None

    Returns:
        Dict of ligand_id:ligand_smiles.
    """
    # Read all ligand data once 
    first_species = next(AF3_INPUT_DIR.iterdir())
    first_protein = next(first_species.iterdir())

    #ligand_id : smiles
    ligand_data: dict[str:str] = {}

    for ligand_file in first_protein.iterdir():
        data = json.load(open(ligand_file))
        ligand_info = data['sequences'][1]['ligand']
        ligand_data[ligand_file.stem] = ligand_info['smiles']

    print(f"Loaded {len(ligand_data)} ligands once")
    # pprint(ligand_data)

    excluded_ligands = {'514', '2169'} #had no complexes for c3p0 (which is what were comparing against)
    ligand_data = {k: v for k, v in ligand_data.items() 
                        if k not in excluded_ligands}

    print(f"Processing {len(ligand_data)} ligands (excluded {len(excluded_ligands)})")
    return ligand_data