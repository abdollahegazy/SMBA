
import re
import csv
import json
from pathlib import Path
from tqdm import tqdm

DATA_DIR = Path("../data")
OUTPUT_CSV = DATA_DIR / "docking_results.csv"


def load_pairs(data_dir: Path = DATA_DIR) -> list[tuple[str, str]]:
    pairs = []
    for protein_dir in sorted(data_dir.iterdir()):
        if not protein_dir.is_dir():
            continue
        ligands_dir = protein_dir / "ligands_prepared"
        if not ligands_dir.exists():
            continue
        for ligand_file in sorted(ligands_dir.iterdir()):
            ligand_id = ligand_file.stem
            pairs.append((protein_dir.name, ligand_id))
    return pairs


def parse_best_affinity(log_file: Path) -> float | None:
    if not log_file.exists() or log_file.stat().st_size == 0:
        return None

    best = None
    in_table = False
    with open(log_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("-----+"):
                in_table = True
                continue
            if in_table:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        affinity = float(parts[1])
                        if best is None or affinity < best:
                            best = affinity
                    except ValueError:
                        in_table = False
                else:
                    in_table = False
    return best


def parse_boltz_affinity(
    protein_id: str, ligand_id: str, data_dir: Path = DATA_DIR
) -> float | None:
    json_path = (
        data_dir
        / protein_id
        / "boltz"
        / ligand_id
        / "boltz_output"
        / "boltz_output"
        / f"boltz_results_{protein_id}_{ligand_id}"
        / "predictions" 
        / f"{protein_id}_{ligand_id}"
        / f"affinity_{protein_id}_{ligand_id}.json"
    )

    json_path1 = (
        data_dir
        / protein_id
        / "boltz"
        / ligand_id
        / "boltz_output"
        / f"boltz_results_{protein_id}_{ligand_id}"
        / "predictions" 
        / f"{protein_id}_{ligand_id}"
        / f"affinity_{protein_id}_{ligand_id}.json"
    )
    if not json_path.exists():
        # print(f"Boltz affinity JSON not found for {protein_id} + {ligand_id}")
        # print(json_path)
        if not json_path1.exists():
            # print(json_path1)
            return None
        else:
            json_path = json_path1
        # return None
    try:
        with open(json_path) as f:
            data = json.load(f)
        y = float(data["affinity_pred_value"])
        return -1*(6 - y) * 1.364
    except (KeyError, ValueError, json.JSONDecodeError):
        return None


def main():
    pairs = load_pairs()
    print(pairs)
    results = []
    for protein_id, ligand_id in tqdm(pairs, desc="Parsing"):
        # --- Vina ---
        docking_dir = DATA_DIR / protein_id / "docking" / ligand_id

        best_affinity = None
        best_conf = None
        best_pocket = None

        if docking_dir.exists():
            for log_file in docking_dir.glob("out_c*_p*.log"):
                match = re.match(r"out_c(\d+)_p(\d+)\.log", log_file.name)
                if not match:
                    continue
                conf_idx = int(match.group(1))
                pocket_idx = int(match.group(2))
                affinity = parse_best_affinity(log_file)
                if affinity is not None and (
                    best_affinity is None or affinity < best_affinity
                ):
                    best_affinity = affinity
                    best_conf = conf_idx
                    best_pocket = pocket_idx

        # --- Boltz ---
        boltz_affinity = parse_boltz_affinity(protein_id, ligand_id, DATA_DIR)

        if best_affinity is not None or boltz_affinity is not None:
            results.append(
                {
                    "protein": protein_id,
                    "ligand": ligand_id,
                    "vina_affinity": best_affinity,
                    "vina_conformation": best_conf,
                    "vina_pocket": best_pocket,
                    "boltz_affinity": boltz_affinity,
                }
            )

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "protein",
                "ligand",
                "vina_affinity",
                "vina_conformation",
                "vina_pocket",
                "boltz_affinity",
            ],
        )
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {len(results)} rows to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()

# import re
# import csv
# from pathlib import Path
# from tqdm import tqdm

# DATA_DIR = Path("../data")
# OUTPUT_CSV = DATA_DIR / "docking_results.csv"


# def load_pairs(data_dir: Path = DATA_DIR) -> list[tuple[str, str]]:
#     """
#     Discover (protein_id, ligand_id) from directory structure.
#     """

#     pairs = []
#     for protein_dir in sorted(data_dir.iterdir()):
#         if not protein_dir.is_dir():
#             continue
#         ligands_dir = protein_dir / "ligands_prepared"
#         if not ligands_dir.exists():
#             continue

#         for ligand_file in sorted(ligands_dir.iterdir()):
#             ligand_id = ligand_file.stem

#             pairs.append((protein_dir.name, ligand_id))
#     return pairs


# def parse_best_affinity(log_file: Path) -> float | None:
#     if not log_file.exists() or log_file.stat().st_size == 0:
#         return None

#     best = None
#     in_table = False
#     with open(log_file) as f:
#         for line in f:
#             line = line.strip()
#             if line.startswith("-----+"):
#                 in_table = True
#                 continue
#             if in_table:
#                 parts = line.split()
#                 if len(parts) >= 2:
#                     try:
#                         affinity = float(parts[1])
#                         if best is None or affinity < best:
#                             best = affinity
#                     except ValueError:
#                         in_table = False
#                 else:
#                     in_table = False
#     return best


# def main():
#     pairs = load_pairs()
#     print(pairs)
#     results = []
#     for protein_id, ligand_id in tqdm(pairs, desc="Parsing"):
#         docking_dir = DATA_DIR / protein_id / "docking" / ligand_id

#         if not docking_dir.exists():
#             continue

#         best_affinity = None
#         best_conf = None
#         best_pocket = None

#         for log_file in docking_dir.glob("out_c*_p*.log"):
#             match = re.match(r"out_c(\d+)_p(\d+)\.log", log_file.name)
#             if not match:
#                 continue

#             conf_idx = int(match.group(1))
#             pocket_idx = int(match.group(2))

#             affinity = parse_best_affinity(log_file)
#             if affinity is not None and (
#                 best_affinity is None or affinity < best_affinity
#             ):
#                 best_affinity = affinity
#                 best_conf = conf_idx
#                 best_pocket = pocket_idx

#         if best_affinity is not None:
#             results.append(
#                 {
#                     "protein": protein_id,
#                     "ligand": ligand_id,
#                     "best_affinity": best_affinity,
#                     "conformation": best_conf,
#                     "pocket": best_pocket,
#                 }
#             )

#     with open(OUTPUT_CSV, "w", newline="") as f:
#         writer = csv.DictWriter(
#             f,
#             fieldnames=["protein", "ligand", "best_affinity", "conformation", "pocket"],
#         )
#         writer.writeheader()
#         writer.writerows(results)

#     print(f"Wrote {len(results)} rows to {OUTPUT_CSV}")


# if __name__ == "__main__":
#     main()
