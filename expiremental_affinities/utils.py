from typing import Dict, Optional
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

def pdb_to_fasta_by_chain(pdb_path: str) -> Optional[Dict[str, str]]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)

    model = next(structure.get_models())
    out: Dict[str, str] = {}

    for chain in model:
        aa = []
        for res in chain:
            # res.id[0] == " " means standard residue (not HETATM/water)
            if res.id[0] != " ":
                continue
            try:
                aa.append(seq1(res.resname))
            except KeyError:
                # non-standard residue name -> skip or set 'X'
                aa.append("X")
        seq = "".join(aa)
        if seq:
            out[chain.id] = seq

    return out if out else None


# seqs = pdb_to_fasta_by_chain("./PDBbind/P-L/1gfz/1gfz_protein.pdb")
# print(seqs.keys())  # real chain IDs from the PDB, e.g. dict_keys(['A'])

seqs = pdb_to_fasta_by_chain("./PDBbind/P-L/6cnj/6cnj_protein.pdb")
print(seqs.keys())
print(seqs["A"][:80])

print(seqs["B"][:80])
