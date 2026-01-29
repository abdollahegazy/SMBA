"""
this checks if two given smiles are the same 
by getting their mol from smiles, removing stereochemistry (bc this might look diff in smiles)
then converting back to smiles and making sure its the same.

just used for fun. (this is quinine)
"""
from rdkit import Chem

chemspi_mol = Chem.MolFromSmiles("C=C[C@H]1CN2CC[C@H]1C[C@H]2[C@H](O)c1ccnc2ccc(OC)cc12")
rcsb_mol = Chem.MolFromSmiles("O(c4cc1c(nccc1C(O)C2N3CCC(C2)C(/C=C)C3)cc4)C")

Chem.RemoveStereochemistry(chemspi_mol)
Chem.RemoveStereochemistry(rcsb_mol)

chemspi_flat = Chem.MolToSmiles(chemspi_mol)
rcsb_flat = Chem.MolToSmiles(rcsb_mol)

print(chemspi_flat == rcsb_flat) 