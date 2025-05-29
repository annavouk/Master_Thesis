import pandas as pd
import pickle

df = pd.read_csv("compound_seed_nonseed_counts_with_kegg.csv", low_memory=False)

PICKLE_FILE_SEEDs = "input/seeds_binary_per_patric.pckl"
PICKLE_FILE_non_SEEDs = "input/non_seeds_binary_per_patric.pckl"

seed_df = pd.read_pickle(PICKLE_FILE_SEEDs)
non_seed_df = pd.read_pickle(PICKLE_FILE_non_SEEDs)

# Get compound IDs from pickles and combine them to a list
seed_cpds = set(seed_df.columns[(seed_df == 1).any()])
nonseed_cpds = set(non_seed_df.columns[(non_seed_df == 1).any()])
cpds_list = list(seed_cpds.union(nonseed_cpds))

# Check if all ModelSEED IDs are in the 'ModelSEED ID' column of the DataFrame
missing_ids = [id for id in cpds_list if id not in df["ModelSEED ID"].values]

if missing_ids:
    print(f"The following ModelSEED IDs are missing: {missing_ids}")
else:
    print("All ModelSEED IDs are present in the DataFrame.")
