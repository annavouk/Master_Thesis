import pandas as pd
import numpy as np

def seeds_per_genome(df):
	# Count 1s (seeds) per row
	df['Total_Seeds'] = np.sum(df.values == 1, axis=1)
	most_compounds_genome = df['Total_Seeds'].idxmax()
	least_compounds_genome = df['Total_Seeds'].idxmin()
	max_compounds = df.loc[most_compounds_genome, 'Total_Seeds']
	min_compounds = df.loc[least_compounds_genome, 'Total_Seeds']

	return (df, most_compounds_genome, max_compounds, least_compounds_genome, min_compounds)

if __name__ == "__main__":
	df = pd.read_pickle('seeds_binary_per_patric.pckl')

	(result_df, most_compounds_genome, max_count, least_compounds_genome, min_count) = seeds_per_genome(df)

	print(result_df)
	print(f"The genome with the most seed compounds: {most_compounds_genome} (count: {max_count})")
	print(f"The genome with the least seed compounds: {least_compounds_genome} (count: {min_count})")

# Save result_df to a pkl file for another use
result_df.to_pickle("seeds_per_genome_result_df.pkl")

