import pandas as pd

df = pd.read_csv('merged_metadata.csv', low_memory=False)

def get_unique_species(df):
    unique_species = []
    species_set = set()

    for taxonomy in df['gtdb_taxonomy'].dropna():
        parts = taxonomy.split(';')
        for part in parts:
            if part.startswith('s__'):
                species_name = part.split('__')[1]
                if species_name not in species_set:
                    species_set.add(species_name)
                    unique_species.append(species_name)

    return unique_species

unique_species_list = get_unique_species(df)

print(unique_species_list)

