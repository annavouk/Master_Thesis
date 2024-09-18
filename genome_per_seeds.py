import pandas as pd

def genome_per_seeds(df):
	genomes = df.sum()
	return genomes

if __name__ == "__main__":

	df = pd.read_pickle('seeds_binary_per_patric.pckl')

	genome_counts = genome_per_seeds(df)

	new_row = pd.DataFrame(genome_counts).T
	new_row.index = ['Total_Genomes_per_Seed']

	updated_df = pd.concat([df, new_row])

	print(updated_df)
