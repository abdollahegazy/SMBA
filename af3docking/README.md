# AlphaFold3 Docking

This directory runs AlphaFold3 (AF3) on the full SMBA protein set, including associated ligands, so that AF3 performs protein–ligand docking alongside structure prediction.

All SMBA proteins were processed with their ligands, resulting in AF3-predicted complex structures for the entire dataset.

From these results, we selected the 20 complexes with the highest ligand binding-site RMSD (after protein alignment) and simulated them in NAMD for 250 ns each.
