# Boltz-2 Predictions

This directory contains large-scale Boltz-2 predictions for all SMBA protein–ligand complexes (~20,000 proteins × 12 ligands), used to obtain predicted structures and binding affinities.

Affinities are compared against the SMBA pipeline. From this comparison, we selected:
- 20 complexes originating from SMBA with affinity ≤ −10 kcal/mol and the largest differences vs Boltz-2  
- 20 complexes originating from Boltz-2 with affinity ≤ −8 kcal/mol and the largest differences vs SMBA  

For each selected complex, both the SMBA and Boltz-2 structures are simulated.  
This results in **80 total complexes** simulated in NAMD, each for **250 ns**.

Boltz-2 affinity scores are converted according to the official documentation:  
https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md#output
