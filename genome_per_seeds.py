import pandas as pd

def genome_per_seeds(df):
	# Count 1s per column
	genomes = df.sum()
	most_genomes = genomes.idxmax()
	least_genomes = genomes.idxmin()
	max_count = genomes.max()
	min_count = genomes.min()

	# Count 0s per row
	non_seeds = (df == 0).sum()
	most_non_seed_genomes = non_seeds.idxmax()
	least_non_seed_genomes = non_seeds.idxmin()
	max_non_seed_count = non_seeds.max()
	min_non_seed_count = non_seeds.min()

	return (genomes, most_genomes, max_count, least_genomes, min_count, non_seeds, most_non_seed_genomes, max_non_seed_count, least_non_seed_genomes, min_non_seed_count)


if __name__ == "__main__":
	df = pd.read_pickle('seeds_binary_per_patric.pckl')
	(genomes, most_genomes, max_count, least_genomes, min_count, non_seeds, most_non_seed_genomes, max_non_seed_count, least_non_seed_genomes, min_non_seed_count) = genome_per_seeds(df)

	sum_row = pd.DataFrame(genomes).T
	sum_row.index = ['Total_Genomes_per_Seed']
	sum_row_non_seeds = pd.DataFrame(non_seeds).T
	sum_row_non_seeds.index = ['Total_Genomes_per_non_seed']
	updated_df = pd.concat([df, sum_row, sum_row_non_seeds])

	print(updated_df)
	print(f"The seed compound present in the most genomes: {most_genomes} (count: {max_count})")
	print(f"The seed compound present in the least genomes: {least_genomes} (count: {min_count})")
	print(f"The non-seed compound present in the most genomes: {most_non_seed_genomes} (count: {max_non_seed_count})")
	print(f"The non-seed compound present in the least genomes: {least_non_seed_genomes} (count: {min_non_seed_count})")
