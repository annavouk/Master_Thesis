import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('seeds_to_non_seeds_genome_size.csv', low_memory=False)

def get_genomes_by_taxonomy(df, level_or_name):
    taxonomic_levels = {
        'domain': 'd__',
        'phylum': 'p__',
        'class': 'c__',
        'order': 'o__',
        'family': 'f__',
        'genus': 'g__',
        'species': 's__'
    }

    # Filter by the genus or taxonomic level
    if level_or_name.lower() in taxonomic_levels:
        prefix = taxonomic_levels[level_or_name.lower()]
        filtered_df = df[df['gtdb_taxonomy'].str.contains(prefix, na=False)]
    else:
        filtered_df = df[df['gtdb_taxonomy'].str.contains(level_or_name, case=False, na=False)]

    return filtered_df[['patric_id','genome_size', 'gtdb_taxonomy', 'Total_Seeds', 'Total_non_Seeds']].dropna()

# Example: Specify the genus or taxonomic level you want to filter by
level_or_name = 'Rhizobium'

# Get the filtered DataFrame
genomes_df = get_genomes_by_taxonomy(df, level_or_name)

# Print the filtered DataFrame to inspect the data
print(genomes_df)

# Check how many genomes are found
print(f"Number of Rhizobium genomes in the dataset: {len(genomes_df)}")

# Plotting the Total Seeds for each individual Rhizobium genome
plt.figure(figsize=(20, 10))
plt.bar(genomes_df['patric_id'].astype(str), genomes_df['Total_Seeds'], color='green')
plt.xlabel('Patric IDs')
plt.ylabel('Total Seeds')
plt.title(f'Seeds in Each Rhizobium Genome')
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
