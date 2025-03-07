import pandas as pd

df1 = pd.read_csv("seeds_to_non_seeds.csv")

df2 = pd.read_csv("compact_metadata.csv")

merged_df = pd.merge(df1, df2[['patric_id', 'genome_size', 'gtdb_taxonomy']], left_on='PATRIC', right_on='patric_id', how='left')

merged_df.to_csv("seeds_to_non_seeds_genome_size.csv", index=False)
