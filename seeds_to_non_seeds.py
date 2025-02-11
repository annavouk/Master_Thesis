import pandas as pd

seeds_df = pd.read_pickle('seeds_per_genome_result_df.pkl')
non_seeds_df = pd.read_pickle('non_seeds_per_genome_result_df.pkl')

merged_df = pd.merge(seeds_df['Total_Seeds'], non_seeds_df['Total_non_Seeds'], on="PATRIC")

merged_df['Ratio'] = seeds_df['Total_Seeds'] / non_seeds_df['Total_non_Seeds']

print(merged_df)

merged_df.to_pickle('seeds_to_non_seeds_df.pkl')

min = merged_df[['Ratio']].min()
max = merged_df[['Ratio']].max()

print(min)
print(max)
