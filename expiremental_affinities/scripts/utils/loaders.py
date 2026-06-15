import csv
from pathlib import Path

DATA_DIR = Path("../data/predictions")
KNOWN_STRUCTURES_CSV = Path("../data/known_structures.csv")

def load_pairs(csv_path: Path = KNOWN_STRUCTURES_CSV) -> list[tuple[str, str, str]]:
    """
    Read (protein_id, ligand_id, smiles) from known_structures.csv.
    protein_id is measurement_id.
    """
    pairs = []
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            pairs.append((row["measurement_id"], row["ligand_id"].split(".")[0], row["smiles"]))
    return pairs

def load_proteins(data_dir: Path = DATA_DIR) -> list[str]:
    """
    Discover protein IDs from directory structure.
    """
    proteins = []
    for protein_dir in sorted(data_dir.iterdir()):
        if not protein_dir.is_dir():
            continue
        proteins.append(protein_dir.name)
    return proteins