import pandas as pd

# PATRIC IDs of each Category
#df_1 = pd.read_csv('category_1_genomes.csv')
#ids_1 = df_1['PATRIC'].astype(str).tolist()

#df_2 = pd.read_csv('category_2_genomes.csv')
#ids_2 = df_2['PATRIC'].astype(str).tolist()

df_3 = pd.read_csv('category_3_genomes.csv')
ids_3 = df_3['PATRIC'].astype(str).tolist()

# Filter the df
def get_seeds_subset(df, patric_ids):
	existing_ids = df.index.intersection(patric_ids)
	subset = df.loc[existing_ids]
	return subset

if __name__ == "__main__":
	seeds_binary_df = pd.read_pickle('seeds_binary_per_patric.pckl')
	#non_seeds_binary_df = pd.read_pickle('non_seeds_binary_per_patric.pckl')

	subset_df = get_seeds_subset(seeds_binary_df, ids_3)

	print("\nSubset:")
	print(subset_df)
