import pandas as pd
import pickle
import re

ModelSEED_FILE = "compounds.tsv"

PICKLE_FILE_SEEDs = "seeds_binary_per_patric.pckl"
PICKLE_FILE_non_SEEDs = "non_seeds_binary_per_patric.pckl"

seed_df = pd.read_pickle(PICKLE_FILE_SEEDs)
non_seed_df = pd.read_pickle(PICKLE_FILE_non_SEEDs)
db_df = pd.read_csv(ModelSEED_FILE, delimiter='\t', low_memory = False)

# Get compound IDs from pickles and combine them to a list
seed_cpds = set(seed_df.columns[(seed_df == 1).any()])
nonseed_cpds = set(non_seed_df.columns[(non_seed_df == 1).any()])
cpds_list = list(seed_cpds.union(nonseed_cpds))

filtered_db_df = db_df[db_df['id'].isin(cpds_list)]

def extract_kegg_id(alias_string):
    if pd.notna(alias_string):
        match = re.search(r'\|KEGG: ([^|]+)', alias_string)
        if match:
            return [k.strip() for k in match.group(1).split(';') if k.strip()]
    return []

filtered_db_df = filtered_db_df.copy()
filtered_db_df['kegg_ids'] = filtered_db_df['aliases'].apply(extract_kegg_id)
id_to_kegg_map = dict(zip(filtered_db_df['id'], filtered_db_df['kegg_ids']))

df = pd.DataFrame(list(id_to_kegg_map.items()), columns=['ModelSEED Compound ID', 'KEGG ID'])
df.to_csv('modelseed_to_kegg.csv', index=False)
