from MDAnalysis.analysis import align
from MDAnalysis import Universe

def rmsd_align(
    moving_universe: Universe,
    reference_universe: Universe,
    moving_sel: str,
    reference_sel: str,
) -> None:
    """
    Aligns the entire moving pdb to the reference pdb based on the specified atom selections, using RMSD minimization.

    Parameters:
    - moving_universe (Universe): MDAnalysis Universe object for the structure to be aligned (the "moving" structure).
    - reference_universe (Universe): MDAnalysis Universe object for the reference structure to which the moving structure will be aligned.
    - moving_sel (str): MDAnalysis atom selection string specifying the subset of atoms in the moving structure to use for alignment (e.g., "name CA" to select alpha carbons).
    - reference_sel (str): MDAnalysis atom selection string specifying the subset of atoms in the reference structure to use for alignment (e.g., "name CA" to select alpha carbons).
    - output_pdb (Path | None): Optional path to write the aligned structure.
    """

    moving = moving_universe
    reference = reference_universe

    old_rmsd, new_rmsd = align.alignto(
        mobile=moving,
        reference=reference,
        select={"mobile": moving_sel, "reference": reference_sel},
    )

    # print(f"RMSD before alignment: {old_rmsd:.3f} Å"f"\nRMSD after alignment: {new_rmsd:.3f} Å")

def ligand_center(universe: Universe, ligand_sel: str) -> tuple[float, float, float]:
    ligand = universe.select_atoms(ligand_sel)
    return ligand.center_of_geometry()