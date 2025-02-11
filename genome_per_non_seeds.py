import pandas as pd

def genome_per_non_seeds(df):
	# Count 1s per column
	genomes = df.sum()
	most_genomes = genomes.idxmax()
	least_genomes = genomes.idxmin()
	max_count = genomes.max()
	min_count = genomes.min()

	return (genomes, most_genomes, max_count, least_genomes, min_count)


if __name__ == "__main__":
	df = pd.read_pickle('non_seeds_binary_per_patric.pckl')
	(genomes, most_genomes, max_count, least_genomes, min_count) = genome_per_non_seeds(df)

	sum_row = pd.DataFrame(genomes).T
	sum_row.index = ['Total_Genomes_per_non_Seed']
	updated_df = pd.concat([df, sum_row])

	print(updated_df)
	print(f"The non seed compound present in the most genomes: {most_genomes} (count: {max_count})")
	print(f"The non seed compound present in the least genomes: {least_genomes} (count: {min_count})")
