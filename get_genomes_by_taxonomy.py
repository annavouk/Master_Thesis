import pandas as pd

df = pd.read_csv('merged_metadata.csv', low_memory=False)

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

    if level_or_name.lower() in taxonomic_levels:
        prefix = taxonomic_levels[level_or_name.lower()]
        filtered_df = df[df['gtdb_taxonomy'].str.contains(prefix, na=False)]
    else:
        filtered_df = df[df['gtdb_taxonomy'].str.contains(level_or_name, case=False, na=False)]

    genomes = filtered_df[['patric_id', 'assembly_accession']].dropna().values.tolist()
    return genomes

# Example:
level_or_name = 'Firmicutes'
#level_or_name = 'species'
genomes_list = get_genomes_by_taxonomy(df, level_or_name)

print(f"Genomes for {level_or_name}:")
for genome in genomes_list:
    print(genome)
