"""Build best-pose protein-ligand complexes for all three docking sources.

One CGenFF param set per ligand lives in
    scripts/simulations/ligand_topos_params/{ligand}/l{NN}.rtf  (+ .prm)
(generated once from SMILES). For every source we relabel the docked ligand's
heavy atoms with that rtf's atom names (RDKit graph match), keeping the pose
coordinates, so build_system.tcl can read the ligand straight from the complex
with no per-source renaming. The same rtf/prm then serves all three sims.

Sources:
  smba             docking/docking/{lig}                  -> conformation protein
  smba_af3_pocket  docking_af3_pocket/docking_af3_pocket/{lig} -> conformation protein
  boltz            boltz/{lig}/boltz_output/.../*_model_0.pdb  -> boltz-predicted protein

Ligands without a ligand_topos rtf are skipped (existing complexes untouched).
"""

import re
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem

from utils import DATA_DIR

POSE_RE = re.compile(r"out_poses_c(\d+)_p(\d+)\.pdbqt$")
COMBINE_TCL = Path(__file__).resolve().parent / "combine.tcl"

# One param set per ligand: ligand_topos_params/{ligand}/l{NN}.rtf (+ .prm)
LIG_TOPO_BASE = Path(__file__).resolve().parent.parent / "simulations" / "ligand_topos_params"

# Vina docking sources, relative to each protein dir.
VINA_SOURCES = [
    ("docking/docking", "smba"),
    ("docking_af3_pocket/docking_af3_pocket", "smba_af3_pocket"),
]


def ligand_rtf(ligand: str) -> Path | None:
    """The single CGenFF rtf for a ligand, or None if not parameterized yet."""
    rtfs = sorted((LIG_TOPO_BASE / ligand).glob("*.rtf"))
    return rtfs[0] if rtfs else None


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


def pdbqt_to_pdb(pdbqt_text: str, dst: Path):
    """Convert a ligand pdbqt block to a PDB file via obabel."""
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "lig.pdbqt"
        src.write_text(pdbqt_text)
        subprocess.run(["obabel", str(src), "-O", str(dst)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _ref_graph_from_rtf(rtf_path: Path):
    """Build a heavy-atom RDKit graph from a CGenFF rtf (atom names + bonds).

    Returns (mol, names) where names[i] is the rtf atom name of mol atom i, or
    (None, None) on failure. Bonds are single (connectivity only); hydrogens and
    H-bonds are dropped to match the heavy-atom pose.
    """
    names, idx_of, bonds = [], {}, []
    for line in rtf_path.read_text().splitlines():
        tok = line.split()
        if not tok:
            continue
        if tok[0] == "ATOM" and len(tok) >= 2:
            name = tok[1]
            elem = re.match(r"[A-Za-z]+", name).group().capitalize()
            if elem == "H":
                continue
            idx_of[name] = len(names)
            names.append((name, elem))
        elif tok[0] in ("BOND", "DOUBLE", "TRIPLE"):
            atoms = tok[1:]
            for k in range(0, len(atoms) - 1, 2):
                bonds.append((atoms[k], atoms[k + 1]))
    if not names:
        return None, None
    rw = Chem.RWMol()
    for _, elem in names:
        rw.AddAtom(Chem.Atom(elem))
    for a, b in bonds:
        if a in idx_of and b in idx_of:
            rw.AddBond(idx_of[a], idx_of[b], Chem.BondType.SINGLE)
    return rw.GetMol(), [n for n, _ in names]


def _heavy_mol(path: Path):
    """RDKit mol of heavy atoms only, atom order preserved from the pdb.

    Bonds are flattened to single so matching is purely connectivity-based --
    obabel/boltz and the rtf perceive aromaticity/bond orders differently, which
    would otherwise break the graph match between identical molecules.
    """
    mol = Chem.MolFromPDBFile(str(path), removeHs=True, sanitize=False)
    if mol is None:
        return None
    rw = Chem.RWMol(mol)
    for idx in sorted((a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1),
                      reverse=True):
        rw.RemoveAtom(idx)
    m = rw.GetMol()
    for b in m.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    return m


def _heavy_lines(pdb_path: Path) -> list[str]:
    """Heavy-atom ATOM/HETATM lines from a pdb, in file order."""
    out = []
    for line in pdb_path.read_text().splitlines():
        if line[:6] in ("ATOM  ", "HETATM"):
            elem = (line[76:78].strip() or line[12:16].strip()[:1])
            if elem != "H":
                out.append(line)
    return out


def _fmt_name(name: str) -> str:
    """Pad an atom name into PDB columns 13-16 (CGenFF style)."""
    return name[:4] if len(name) >= 4 else f" {name}".ljust(4)


def name_ligand_lines(pose_pdb: Path, rtf_path: Path) -> list[str] | None:
    """Heavy-atom PDB lines for the pose relabeled with the rtf's CGenFF names
    (graph match), keeping pose coordinates, chain L. None on failure."""
    ref, refnames = _ref_graph_from_rtf(rtf_path)
    pose = _heavy_mol(pose_pdb)
    if ref is None or refnames is None or pose is None:
        return None
    if ref.GetNumAtoms() != pose.GetNumAtoms():
        return None
    match = pose.GetSubstructMatch(ref)  # pose atom idx for each ref atom
    if len(match) != ref.GetNumAtoms():
        return None
    name_of = {match[i]: refnames[i] for i in range(len(match))}

    heavy_lines = _heavy_lines(pose_pdb)
    if len(heavy_lines) != pose.GetNumAtoms():
        return None
    out = []
    for i, line in enumerate(heavy_lines):
        nm = name_of.get(i)
        if nm is None:
            return None
        out.append(line[:12] + _fmt_name(nm) + line[16:21] + "L" + line[22:])
    return out


def _write_pdb(lines: list[str], path: Path):
    path.write_text("\n".join(lines) + "\nEND\n")


def create_complex(protein_pdb: Path, ligand_pdb: Path, outpath: Path):
    """Merge protein (chain P) + ligand (chain L) with topotools via VMD."""
    subprocess.run(
        ["vmd2", "-dispdev", "text", "-e", str(COMBINE_TCL),
         "-args", str(protein_pdb), str(ligand_pdb), str(outpath)],
        check=True,
    )


def boltz_pdb_path(protein_dir: Path, ligand: str) -> Path:
    pid = protein_dir.name
    return (protein_dir / "boltz" / ligand / "boltz_output"
            / f"boltz_results_{pid}_{ligand}" / "predictions"
            / f"{pid}_{ligand}" / f"{pid}_{ligand}_model_0.pdb")


def build_vina(protein_dir: Path, ligand_dir: Path, source: str, rtf: Path):
    """smba / smba_af3_pocket: best vina pose + matching conformation protein."""
    ligand = ligand_dir.name
    tag = f"{protein_dir.name}/{source}/{ligand}"

    best = best_pose_for_ligand(ligand_dir)
    if best is None:
        print(f"[skip] {tag}: no poses")
        return
    pose_file, conformation, affinity = best

    protein_pdb = protein_dir / "conformation" / f"protein_conf{conformation}.pdb"
    if not protein_pdb.exists():
        print(f"[skip] {tag}: missing protein {protein_pdb}")
        return

    pose_pdb = ligand_dir / f"{ligand}_pose.pdb"
    named_pdb = ligand_dir / f"{ligand}_named.pdb"
    outpath = ligand_dir / f"{ligand}_best_complex.pdb"

    pdbqt_to_pdb(extract_model1(pose_file), pose_pdb)
    named = name_ligand_lines(pose_pdb, rtf)
    if named is None:
        print(f"[FAIL] {tag}: could not match pose to {rtf.name}")
        pose_pdb.unlink(missing_ok=True)
        return
    _write_pdb(named, named_pdb)
    create_complex(protein_pdb, named_pdb, outpath)
    pose_pdb.unlink(missing_ok=True)
    named_pdb.unlink(missing_ok=True)
    print(f"[ok] {tag}  conf={conformation}  affinity={affinity}  -> {outpath}")


def build_boltz(protein_dir: Path, ligand: str, rtf: Path):
    """boltz: relabel the boltz-predicted complex's ligand to the rtf names."""
    tag = f"{protein_dir.name}/boltz/{ligand}"
    src = boltz_pdb_path(protein_dir, ligand)
    if not src.exists():
        print(f"[skip] {tag}: missing boltz output {src}")
        return

    protein_lines, ligand_lines = [], []
    for line in src.read_text().splitlines():
        if line[:6] in ("ATOM  ", "HETATM"):
            (ligand_lines if line[21] == "L" else protein_lines).append(line)
    if not ligand_lines:
        print(f"[skip] {tag}: no chain-L ligand in boltz output")
        return

    out_dir = protein_dir / "boltz" / ligand
    pose_pdb = out_dir / f"{ligand}_pose.pdb"
    _write_pdb(ligand_lines, pose_pdb)

    named = name_ligand_lines(pose_pdb, rtf)
    pose_pdb.unlink(missing_ok=True)
    if named is None:
        print(f"[FAIL] {tag}: could not match boltz ligand to {rtf.name}")
        return

    outpath = out_dir / f"{ligand}_best_complex.pdb"
    _write_pdb(protein_lines + named, outpath)
    print(f"[ok] {tag}  -> {outpath}")


def main():
    for protein_dir in sorted(DATA_DIR.iterdir()):
        if not protein_dir.is_dir():
            continue

        # vina sources
        for source_rel, source in VINA_SOURCES:
            docking_root = protein_dir / source_rel
            if not docking_root.is_dir():
                continue
            for ligand_dir in sorted(docking_root.iterdir()):
                if not ligand_dir.is_dir():
                    continue
                rtf = ligand_rtf(ligand_dir.name)
                if rtf is None:
                    print(f"[skip] {protein_dir.name}/{source}/{ligand_dir.name}: "
                          f"no ligand_topos rtf")
                    continue
                build_vina(protein_dir, ligand_dir, source, rtf)

        # boltz source
        boltz_root = protein_dir / "boltz"
        if boltz_root.is_dir():
            for ligand_dir in sorted(boltz_root.iterdir()):
                if not ligand_dir.is_dir():
                    continue
                rtf = ligand_rtf(ligand_dir.name)
                if rtf is None:
                    print(f"[skip] {protein_dir.name}/boltz/{ligand_dir.name}: "
                          f"no ligand_topos rtf")
                    continue
                build_boltz(protein_dir, ligand_dir.name, rtf)


if __name__ == "__main__":
    main()
