from pathlib import Path

DATA_DIR = Path("../data")

def load_pairs(data_dir: Path = DATA_DIR) -> list[tuple[str, str, str]]:
    """
    Discover (protein_id, ligand_id, smiles_or_ccd) from directory structure.
    """

    pairs = []
    for protein_dir in sorted(data_dir.iterdir()):
        if not protein_dir.is_dir():
            continue
        ligands_dir = protein_dir / "ligands"
        if not ligands_dir.exists():
            continue
        
        for lig_file in sorted(ligands_dir.iterdir()):
            ligand_id = lig_file.stem
            if lig_file.suffix == ".smiles":
                value = lig_file.read_text().strip()
            elif lig_file.suffix == ".ccd":
                value = lig_file.read_text().strip()
            else:
                continue
            pairs.append((protein_dir.name, ligand_id, value))
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