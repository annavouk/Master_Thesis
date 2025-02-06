import pandas as pd

df = pd.read_pickle('seeds_binary_per_patric.pckl')
compounds_ids = df.iloc[0, 1:]
print(compounds_ids)

seed_ids_df = pd.read_csv('seedId_keggId_module.tsv', sep='\t')
print(seed_ids_df.head())
