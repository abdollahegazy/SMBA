from rcsbapi.data import DataQuery as Query
from rcsbapi.data import DataQuery
from typing import List,Optional

import re
from typing import List
import pandas as pd
import requests

def get_complex_fasta(rcsb_code: str) -> Optional[str]:
    """Returns the FASTA if exactly one protein chain, else None."""
    url = f"https://www.rcsb.org/fasta/entry/{rcsb_code}/display"
    response = requests.get(url)

    if response.status_code != 200:
        return None    
    fasta = response.text
    num_chains = fasta.count(">")

    if num_chains == 1:
        return fasta
    return None

def find_pdbbind_entries(ligand_rcsb_codes: List[str], pdbbind_id_file: str = "pdbbing_r12020_index.lst") -> pd.DataFrame:
    hits = []
    pattern = re.compile(r'\((' + '|'.join(re.escape(code) for code in ligand_rcsb_codes) + r')\)')
    with open(pdbbind_id_file) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                ligand = m.group(1)  # Capture the matched ligand code
                # Split on // first
                parts = line.strip().split("//")
                left = parts[0].strip()  # "2tpi  2.10  1982  Kd=49uM"
                right = parts[1].strip() if len(parts) > 1 else ""  # "2tpi.pdf (2-mer)"
                
                # Split left side on 2+ spaces
                cols = re.split(r'\s{2,}', left)
                cols.append(right)
                cols.append(ligand)  # Add ligand as last column
                hits.append(cols)
    
    structures = pd.DataFrame(hits, columns=["PDB ID", "resolution", "year", "binding data", "comments", "ligand"])
    return structures

structures = find_pdbbind_entries(["QI9","NCT","CFF"])
# Filter to only single-protein structures
structures["fasta"] = structures["PDB ID"].apply(get_complex_fasta)
useable = structures[structures["fasta"].notna()]

print(f"Useable: {len(useable)} / {len(structures)}")
print(useable)

# structures = find_pdbbind_entries(["QI9","NCT","CFF"])
# useable = 0
# for structure in structures.iterrows():
#     structure_rcsb_code = structure[1]["PDB ID"]
#     if not get_complex_fasta(structure_rcsb_code):
        
#     print("\n\n\n")
# print(useable)
