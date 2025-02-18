import pandas as pd

compounds_df = pd.read_csv("unique_compounds.csv")

seed_ids_df = pd.read_csv('seedId_keggId_module.tsv', sep='\t')
print(seed_ids_df.head())
