import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the seeds binary matrix: rows = genomes, columns = compounds (1 = seed, 0 = not)
seeds_df = pd.read_pickle("input/seeds_binary_per_patric.pckl")

# For each compound, count in how many genomes it is a seed
seed_freq = seeds_df.sum(axis=0).sort_values(ascending=False)

# Summary: Top-N essential compounds (most frequent seeds)
top_n = 20
top_seeds = seed_freq.head(top_n)

# 1. Barplot: Top-N most frequent seeds
plt.figure(figsize=(12, 5))
ax = sns.barplot(x=top_seeds.index, y=top_seeds.values, color="teal")
plt.xticks(rotation=45, ha="right")
plt.xlabel("Compound (ModelSEED ID)")
plt.ylabel("Number of Genomes with a specific Compound as Seed")
plt.title(f"Top {top_n} Most Essential Metabolites (Seed Frequency)")
for i, v in enumerate(top_seeds.values):
    ax.text(i, v+0.01*max(top_seeds.values), str(v), ha='center', va='bottom', fontsize=9)
plt.tight_layout()
plt.show()

# 2. Histogram: Distribution of seed frequencies
plt.figure(figsize=(7, 4))
plt.hist(seed_freq.values, bins=30, color="slateblue", edgecolor="black", log=True)
plt.xlabel("Number of Genomes with a Specific Compound as Seed")
plt.ylabel("Number of Compounds (log scale)")
plt.title("Distribution of Seed Node Frequencies Across All Compounds")
plt.tight_layout()
plt.show()

print(seed_freq.describe())
