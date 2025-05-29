import pandas as pd
from metabolic_potential_visualization_1 import split_and_clean_taxonomy

# Load the data
df = pd.read_csv("../analysis/metabolic_potential_summary.csv", low_memory=False)

taxonomy_df = split_and_clean_taxonomy(df, 'gtdb_taxonomy')

df = pd.concat([df, taxonomy_df], axis=1)

num_phyla = df['phylum'].nunique()
print(f"Total number of different phyla: {num_phyla}")
