"""
this looks at the PDBBind list of data
and searche for ligands specified in the dict below
for all complexes that have it

then it converts that pdb to fasta and gets ligand smiles
"""

import sys
sys.path.append("/home/abdolla/anaconda3/lib/python3.12/site-packages")

from typing import List
import re
import pandas as pd
import os
import subprocess
from dotenv import load_dotenv
from chemspipy import ChemSpider
from vmd import molecule, evaltcl,atomsel

from pprint import pprint

rcsb2chemspi = {
    "NCT":917,
    "CFF": 2424,
    "QI9": 84989
}




def pdb_to_fasta_via_obabel(pdb_path: str) -> str:
    molid = molecule.load("pdb",pdb_path )
    evaltcl(f"set prot [atomselect {molid} \"protein\"]")
    unique_chains = evaltcl("lsort -unique [$prot get chain]"); n_chains = len(unique_chains.split())
    absolute_pdb_dir = os.path.dirname(os.path.abspath(pdb_path))
    fasta = {}

    for chain in unique_chains.split():
        evaltcl(f'set sel [atomselect {molid} "chain {chain}"]')
        #vmd needs absolute path otherwise did it in cwd
        chain_pdb = f"{absolute_pdb_dir}/chain_{chain}.pdb"

        evaltcl(f'$sel writepdb "{absolute_pdb_dir}/chain_{chain}.pdb"')
        evaltcl("$sel delete")
        res = subprocess.run(
            ["obabel", chain_pdb, "-ofasta"],
            capture_output=True,
            text=True
        )

        fasta[chain] = "".join(res.stdout.split("\n")[1:])
        os.remove(chain_pdb)

    evaltcl("$prot delete")
    molecule.delete(molid)

    return fasta, n_chains

def get_complex_fasta(rcsb_code: str, base_dir: str = "./PDBbind/P-L") -> str:
    rcsb_code = rcsb_code.lower()
    pdb_path = os.path.join(base_dir, rcsb_code, f"{rcsb_code}_protein.pdb")
    fasta, n_chains = pdb_to_fasta_via_obabel(pdb_path)
    return pd.Series({"fasta": fasta, "n_chains": n_chains})


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
    
    structures = pd.DataFrame(hits, columns=["RCSB ID", "resolution", "year", "binding data", "comments", "ligand rcsb"])
    return structures

def get_ligand_smiles(structures: pd.DataFrame, cs: ChemSpider, rcsb2chemspi: dict):
    smiles_cache = {} 
    chem_ids = []
    smiles_list = []

    for _,row in structures.iterrows():
        ligand = row["ligand rcsb"]
        chem_id = rcsb2chemspi.get(ligand)
        chem_ids.append(chem_id)

        if ligand in smiles_cache:
            smiles_list.append(smiles_cache[ligand])
            continue

        details = cs.get_details(chem_id)
        smi = details.get("smiles")

        smiles_cache[ligand] = smi
        smiles_list.append(smi)

    structures["chemspider_id"] = chem_ids
    structures["ligand_smiles"] = smiles_list
    return structures


if __name__ == "__main__":
    load_dotenv()
    CHEMSPI_KEY = os.environ.get("CHEMSPI_KEY")
    CS = ChemSpider(CHEMSPI_KEY)
    structures = find_pdbbind_entries(list(rcsb2chemspi.keys()))
    # get_complex_fasta("5o87")
    # exit()
    structures[["fasta", "n_chains"]] = structures["RCSB ID"].apply(get_complex_fasta)
    structures = get_ligand_smiles(structures,CS,rcsb2chemspi)
    print(structures.drop(columns=["fasta"]))
    structures.to_csv("known_structures.csv",index=False)
    exit()