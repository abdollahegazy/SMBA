import os
import json
from chemspipy import ChemSpider
from pprint import pprint
from tqdm import tqdm  
from copy import deepcopy

cs = ChemSpider('TzfV80skfR1yCAe7oY4y06Q2s5Kzlff4a1NTNVq2')

LIGANDS = [6309, 780, 4444606, 6085, 2169, 23208, 4444418, 6223, 3698, 395716,
           30999, 21105998, 388386, 13876103, 917, 2424, 4450907, 7930, 10194105, 84989, 108426, 393012,
           4518347, 514]

LIGAND_SMILES = {}
LIGANDS = [84989]
from pprint import pprint
print("Fetching SMILES strings...")
# LIGANDS = [7930]
for ligand in tqdm(LIGANDS, desc="Ligands"):
    # pprint(cs.get_details(ligand))

    # pprint(cs.get_datasources())
    LIGAND_SMILES[ligand] = cs.get_details(ligand)['smiles']
    break

print("=" * 5, "LIGAND INFORMATION","="*5)
pprint(LIGAND_SMILES)
print("=" * 10)