"""
"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ast

df = pd.read_csv("/home/annavouk/master_thesis/final_repo/scripts/analysis/compounds_seed_nonseed_with_kegg.csv", low_memory = False)

# Total number of genomes
total_genomes = 33755

# Total number of seed compounds
total_seeds = 452

# Total number of non-seed compounds
total_non_seeds = 545

"""
This script calculates the percentage of genomes in which each metabolite appears as a seed or non-seed.
It then visualizes the distribution of these coverage percentages to compare how widely seed and
non-seed metabolites are shared across genomes.
"""
# Calculate coverage percentage
df['seed_coverage'] = df['seed_genome_count'] / total_genomes * 100
df['non_seed_coverage'] = df['non_seed_genome_count'] / total_genomes * 100

# Set up the plot
plt.figure(figsize=(10, 6))

# Plot the distributions
sns.histplot(df['seed_coverage'], bins=50, kde=False, color='green', label='Seed', alpha=0.6)
sns.histplot(df['non_seed_coverage'], bins=50, kde=False, color='orange', label='Non-Seed', alpha=0.6)

plt.xlabel('Percentage of genomes in which a metabolite is a seed or non-seed')
plt.ylabel('Number of metabolites')
plt.yscale('log')
plt.title('Distribution of seed and non-seed compound coverage')
plt.legend()
plt.tight_layout()
plt.show()
