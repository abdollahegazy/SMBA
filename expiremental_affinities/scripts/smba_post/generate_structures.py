import re
import subprocess
import tempfile
from pathlib import Path

from utils import DATA_DIR

POSE_RE = re.compile(r"out_poses_c(\d+)_p(\d+)\.pdbqt$")


def model1_affinity(pose_file: Path) -> float | None:
    """Return the affinity of MODEL 1 (the best pose) in a vina out_poses file."""
    in_model1 = False
    for line in pose_file.read_text().splitlines():
        if line.startswith("MODEL"):
            in_model1 = int(line.split()[1]) == 1
        elif in_model1 and line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def extract_model1(pose_file: Path) -> str:
    """Return the text block of MODEL 1 only (most negative affinity)."""
    lines = []
    in_model1 = False
    for line in pose_file.read_text().splitlines():
        if line.startswith("MODEL"):
            in_model1 = int(line.split()[1]) == 1
        if in_model1:
            lines.append(line)
        if in_model1 and line.startswith("ENDMDL"):
            break
    return "\n".join(lines) + "\n"


COMBINE_TCL = Path(__file__).resolve().parent / "combine.tcl"


def pdbqt_to_pdb(pdbqt_text: str, dst: Path):
    """Convert a ligand pdbqt block to a PDB file via obabel."""
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "lig.pdbqt"
        src.write_text(pdbqt_text)
        subprocess.run(["obabel", str(src), "-O", str(dst)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def create_complex(protein_pdb: Path, ligand_pdb: Path, outpath: Path):
    """Merge protein (chain P) + ligand (chain L) with topotools via VMD."""
    subprocess.run(
        ["vmd2", "-dispdev", "text", "-e", str(COMBINE_TCL),
         "-args", str(protein_pdb), str(ligand_pdb), str(outpath)],
        check=True,
    )


def best_pose_for_ligand(ligand_dir: Path) -> tuple[Path, int, float] | None:
    """Find the out_poses file whose MODEL 1 has the most negative affinity."""
    best = None  # (pose_file, conformation, affinity)
    for pose_file in ligand_dir.glob("out_poses_c*_p*.pdbqt"):
        m = POSE_RE.search(pose_file.name)
        if not m:
            continue
        affinity = model1_affinity(pose_file)
        if affinity is None:
            continue
        conformation = int(m.group(1))
        if best is None or affinity < best[2]:
            best = (pose_file, conformation, affinity)
    return best


def main():
    for protein_dir in sorted(DATA_DIR.iterdir()):
        docking_root = protein_dir / "docking" / "docking"
        if not docking_root.is_dir():
            continue

        for ligand_dir in sorted(docking_root.iterdir()):
            if not ligand_dir.is_dir():
                continue
            ligand = ligand_dir.name

            best = best_pose_for_ligand(ligand_dir)
            if best is None:
                print(f"[skip] no poses for {protein_dir.name}/{ligand}")
                continue
            pose_file, conformation, affinity = best

            outpath = ligand_dir / f"{ligand}_best_complex.pdb"
            if outpath.exists() and len(outpath.read_text().splitlines()) > 10:
                continue

            protein_pdb = protein_dir / "conformation" / f"protein_conf{conformation}.pdb"
            if not protein_pdb.exists():
                print(f"[skip] missing protein {protein_pdb}")
                continue

            ligand_pdb = ligand_dir / f"{ligand}_best_pose.pdb"
            pdbqt_to_pdb(extract_model1(pose_file), ligand_pdb)
            create_complex(protein_pdb, ligand_pdb, outpath)
            ligand_pdb.unlink()
            print(f"[ok] {protein_dir.name}/{ligand}  conf={conformation}  "
                  f"affinity={affinity}  -> {outpath}")


if __name__ == "__main__":
    main()
