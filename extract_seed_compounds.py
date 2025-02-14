import pandas as pd

def extract_seed_compounds(file_path):
	df = pd.read_pickle(file_path)

	genome_compounds = {}

	for genome_id, row in df.iterrows():
		compounds = row[row == 1].index.tolist()

		if compounds:
			genome_compounds[genome_id] = compounds

	return genome_compounds

file_path = "seeds_binary_per_patric.pckl"
output_file = "seed_compounds_for_each_patric_id.csv"

result = extract_seed_compounds(file_path)

result_df = pd.DataFrame(list(result.items()), columns=['Genome ID', 'Seed Compounds'])
result_df.to_csv(output_file, index=False)
