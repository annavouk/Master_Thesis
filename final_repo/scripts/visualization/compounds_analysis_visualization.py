import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("/home/annavouk/master_thesis/final_repo/scripts/analysis/compound_seed_nonseed_counts_with_kegg.csv", low_memory=False)

# Group by KEGG Module and sum seed_genome_count
module_seed_counts = df.groupby("KEGG Module")["seed_genome_count"].sum().reset_index()

# Sort by total seed count, descending
module_seed_counts = module_seed_counts.sort_values(by="seed_genome_count", ascending=False)

#print(module_seed_counts.head(10))
#print(module_seed_counts.tail(10))

# Group by KEGG Module and sum non_seed_genome_count
module_non_seed_counts = df.groupby("KEGG Module")["non_seed_genome_count"].sum().reset_index()

# Sort by total seed count, descending
module_non_seed_counts = module_non_seed_counts.sort_values(by="non_seed_genome_count", ascending=False)

#print(module_non_seed_counts.head(10))
#print(module_non_seed_counts.tail(10))

# Select top 10 modules for seeds and non-seeds
top_seed = module_seed_counts.head(10)
top_non_seed = module_non_seed_counts.head(10)

# Plot top 10 KEGG Modules by seed genome counts
plt.figure(figsize=(12, 6))
sns.barplot(data=top_seed, y="KEGG Module", x="seed_genome_count", color="steelblue")
plt.title("Top 10 KEGG Modules by Seed Genome Count")
plt.xlabel("Seed Genome Count")
plt.ylabel("KEGG Module")
plt.tight_layout()
plt.show()

# Plot top 10 KEGG Modules by non-seed genome counts
plt.figure(figsize=(12, 6))
sns.barplot(data=top_non_seed, y="KEGG Module", x="non_seed_genome_count", color="seagreen")
plt.title("Top 10 KEGG Modules by Non-Seed Genome Count")
plt.xlabel("Non-Seed Genome Count")
plt.ylabel("KEGG Module")
plt.tight_layout()
plt.show()
