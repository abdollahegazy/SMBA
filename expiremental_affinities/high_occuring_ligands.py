#!/usr/bin/env python3
"""
Analyze ligand occurrences in PDBbind data file.
Extracts ligand codes from parentheses and counts their frequency.
"""

import re
from collections import Counter

def analyze_ligands(filename):
    ligand_pattern = re.compile(r'\(([A-Za-z0-9-]+)\)')
    ligands = []
    
    # Affinity type counters
    affinity_types = Counter()
    affinity_pattern = re.compile(r'(Kd|Ki|IC50|EC50|Ka|Km)[=<>]', re.IGNORECASE)
    
    with open(filename, 'r') as f:
        for line in f:
            # Skip comment lines and empty lines
            if line.startswith('#') or not line.strip():
                continue
            
            # Count affinity types
            aff_match = affinity_pattern.search(line)
            if aff_match:
                aff_type = aff_match.group(1).lower()
                # Normalize case
                if aff_type == 'kd':
                    affinity_types['Kd'] += 1
                elif aff_type == 'ki':
                    affinity_types['Ki'] += 1
                elif aff_type == 'ic50':
                    affinity_types['IC50'] += 1
                elif aff_type == 'ec50':
                    affinity_types['EC50'] += 1
                elif aff_type == 'ka':
                    affinity_types['Ka'] += 1
                elif aff_type == 'km':
                    affinity_types['Km'] += 1
                else:
                    affinity_types[aff_type] += 1
            
            # Find all ligand codes in parentheses
            matches = ligand_pattern.findall(line)
            # Only keep 3-letter codes
            for m in matches:
                if len(m) == 3:
                    ligands.append(m)
    
    # Count occurrences
    counts = Counter(ligands)
    
    # Print affinity type summary
    print("="*50)
    print("Affinity type breakdown:")
    print("="*50)
    total_aff = sum(affinity_types.values())
    for aff_type, count in affinity_types.most_common():
        pct = 100 * count / total_aff if total_aff > 0 else 0
        print(f"{aff_type:>10}: {count:>6} ({pct:5.1f}%)")
    print(f"{'TOTAL':>10}: {total_aff:>6}")
    
    # Print results
    print(f"\nTotal ligand entries: {len(ligands)}")
    print(f"Unique ligands: {len(counts)}")
    print("\n" + "="*50)
    print("Top 30 most common ligands:")
    print("="*50)
    
    for ligand, count in counts.most_common(30):
        print(f"{ligand:>10}: {count:>5} occurrences")

    return counts, affinity_types

if __name__ == "__main__":
    counts, affinity_types = analyze_ligands("pdbbing_r12020_index.lst")