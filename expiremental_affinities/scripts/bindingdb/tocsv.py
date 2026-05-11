import pandas as pd
from pathlib import Path
import ast

HITS = Path("../../data/bindingdb/hits")
OUT_CSV = Path("../../data/known_structures_test.csv")





def process_hit(filename:str):
    data = pd.read_csv(filename)
    affinities = data[['ki_nM', 'ic50_nM', 'kd_nM', 'ec50_nM']]

    # find number of affinitiy values that contain > or < 
    # bad_affinities = affinities.map(lambda x: isinstance(x, str) and ('>' in x or '<' in x)).sum().sum()
    # keep only good affinity values
    affinities = affinities.map(lambda x: x if not (isinstance(x, str) and ('>' in x or '<' in x)) else None)
    affintiy_type = affinities.apply(lambda row: row.first_valid_index(), axis=1)
    affinity_value = affinities.apply(lambda row: row.first_valid_index() and row[row.first_valid_index()], axis=1)

    uid = data.measurement_id
    smiles = data.smiles
    sequence = data.sequences.map(lambda x: ast.literal_eval(x)[0].strip().lstrip())
    ligand_name,ligand_id = filename.name.split("_")


    meta_data = data.drop(["ki_nM", "ic50_nM", "kd_nM", "ec50_nM", "sequences",'smiles','measurement_id'], axis=1)
    meta_data['affinity_value'] = affinity_value
    meta_data['affinity_type'] = affintiy_type

    valid_df = pd.concat([uid, sequence, smiles], axis=1)
    valid_df['ligand_name'], valid_df['ligand_id'] = ligand_name, ligand_id
    
    valid_df['meta_data'] = meta_data.to_dict(orient='records')

    return valid_df


if __name__ == "__main__":
    results = []
    for hit_file in HITS.glob("*.csv"):
        result = process_hit(hit_file)
        results.append(result)
    
    final_df = pd.concat(results, ignore_index=True)
    final_df.to_csv(OUT_CSV, index=False)

