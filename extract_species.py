import pandas as pd
import matplotlib.pyplot as plt

#df = pd.read_csv("low_outliers_gtdb_taxonomy.csv", low_memory = False)
df = pd.read_csv("high_outliers_gtdb_taxonomy.csv", low_memory = False)

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

species_counts = species_df['Species'].value_counts()

total_species = species_df['Species'].nunique()
print(f"Total number of unique species: {total_species}")
