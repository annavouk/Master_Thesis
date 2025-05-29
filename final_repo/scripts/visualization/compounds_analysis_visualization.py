import pandas as pd

df = pd.read_csv("~/master_thesis/final_repo/compound_seed_nonseed_counts_with_kegg.csv", low_memory = False)

# Group by KEGG Module and sum seed_genome_count
module_seed_counts = df.groupby("KEGG Module")["seed_genome_count"].sum().reset_index()

# Sort by total seed count, descending
module_seed_counts = module_seed_counts.sort_values(by="seed_genome_count", ascending=False)

print(module_seed_counts.head())
print(module_seed_counts.tail())

# Group by KEGG Module and sum non_seed_genome_count
module_non_seed_counts = df.groupby("KEGG Module")["non_seed_genome_count"].sum().reset_index()

# Sort by total seed count, descending
module_non_seed_counts = module_non_seed_counts.sort_values(by="non_seed_genome_count", ascending=False)

print(module_non_seed_counts.head())
print(module_non_seed_counts.tail())
