import pandas as pd
import numpy as np

def seeds_per_genome(df):
	# Count 1s (seeds) per row
	df['Total_Seeds'] = np.sum(df.values == 1, axis=1)
	most_compounds_genome = df['Total_Seeds'].idxmax()
	least_compounds_genome = df['Total_Seeds'].idxmin()
	max_compounds = df.loc[most_compounds_genome, 'Total_Seeds']
	min_compounds = df.loc[least_compounds_genome, 'Total_Seeds']

	# Count 0s (non-seeds) per row
	df['Total_non_Seeds'] =  np.sum(df.values == 0, axis=1)
	most_non_seeds_genome = df['Total_non_Seeds'].idxmax()
	least_non_seeds_genome = df['Total_non_Seeds'].idxmin()
	max_non_seeds = df.loc[most_non_seeds_genome, 'Total_non_Seeds']
	min_non_seeds = df.loc[least_non_seeds_genome, 'Total_non_Seeds']

	# Calculate the ratio Seeds/non_Seeds
	df['Seed_to_NonSeed_Ratio'] = df['Total_Seeds'] / df['Total_non_Seeds']

	return (df, most_compounds_genome, max_compounds, least_compounds_genome, min_compounds,
            most_non_seeds_genome, max_non_seeds, least_non_seeds_genome, min_non_seeds)

if __name__ == "__main__":
	df = pd.read_pickle('seeds_binary_per_patric.pckl')

	(result_df, most_genome, max_count, least_genome, min_count, most_non_seeds_genome, max_non_seeds, least_non_seeds_genome, min_non_seeds) = seeds_per_genome(df)

	print(result_df)
	print(f"The genome with the most compounds: {most_genome} (count: {max_count})")
	print(f"The genome with the least compounds: {least_genome} (count: {min_count})")
	print(f"The genome with the most non seeds compounds: {most_non_seeds_genome} (count: {max_non_seeds})")
	print(f"The genome with the least non seeds compounds: {least_non_seeds_genome} (count: {min_non_seeds})")
	print(f"The range of the ratio Seeds/non_Seeds is: {min_count / max_non_seeds:.2f} to {max_count / min_non_seeds:.2f}")

# Save result_df to a pkl file for another use
result_df.to_pickle("seeds_per_genome_result_df.pkl")

