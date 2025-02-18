import pandas as pd

df = pd.read_csv('merged_metadata.csv', low_memory=False)

def extract_species(df):
	species = {
	part.split('__')[1]
	for taxonomy in df['gtdb_taxonomy'].dropna()
	for part in taxonomy.split(';')
	if part.startswith('s__')
	}

	return sorted(species)

species_list = extract_species(df)

species_df = pd.DataFrame(species_list, columns=['Species'])

print(species_df.head())
print(species_df.tail())

