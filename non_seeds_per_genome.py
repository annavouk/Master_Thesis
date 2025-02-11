import pandas as pd
import numpy as np

def non_seeds_per_genome(df):
	# Count 1s (non_seeds) per row
	df['Total_non_Seeds'] = np.sum(df.values == 1, axis=1)
	most_non_seed_genome = df['Total_non_Seeds'].idxmax()
	least_non_seed_genome = df['Total_non_Seeds'].idxmin()
	max_non_seeds = df.loc[most_non_seed_genome, 'Total_non_Seeds']
	min_non_seeds = df.loc[least_non_seed_genome, 'Total_non_Seeds']

	return (df, most_non_seed_genome, max_non_seeds, least_non_seed_genome, min_non_seeds)

if __name__ == "__main__":
	df = pd.read_pickle('non_seeds_binary_per_patric.pckl')

	(result_df, most_non_seed_genome, max_non_seeds, least_non_seed_genome, min_non_seeds) = non_seeds_per_genome(df)

	print(result_df)
	print(f"The genome with the most non seeds: {most_non_seed_genome} (count: {max_non_seeds})")
	print(f"The genome with the least non seeds: {least_non_seed_genome} (count: {min_non_seeds})")

# Save result_df to a pkl file for another use
result_df.to_pickle("non_seeds_per_genome_result_df.pkl")

