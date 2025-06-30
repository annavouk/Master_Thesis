"""

"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
#from utils import load_data
from config import OUTPUT_DIR

# Histogram: Distribution of KEGG pathways per seed compound
df_pathways_per_seed = pd.read_csv(OUTPUT_DIR / "counts_pathways_per_seed.csv")
plt.figure(figsize=(6,4))
df_pathways_per_seed['num_kegg_pathways'].hist(bins=20, color='teal', edgecolor='black')
plt.xlabel('Number of KEGG Pathways per Seed')
plt.ylabel('Number of Seed Compounds')
plt.title('Distribution of KEGG Pathways per Seed')
plt.tight_layout()
plt.show()

# Barplot: Seeds per KEGG pathway (top 20)
df_seeds_per_pathway = pd.read_csv(OUTPUT_DIR / "counts_seeds_per_pathway.csv")
df_seeds_per_pathway.columns = ['KEGG_pathway', 'num_seeds']
topN = 20
plt.figure(figsize=(10,5))
sns.barplot(data=df_seeds_per_pathway.head(topN), x='KEGG_pathway', y='num_seeds', color='steelblue')
plt.xticks(rotation=45, ha='right')
plt.xlabel('KEGG Pathway')
plt.ylabel('Number of Seed Compounds')
plt.title(f'Top {topN} KEGG Pathways by Number of Seed Compounds')
plt.tight_layout()
plt.show()

# Barplot: Seeds per KEGG module (top 20)
df_seeds_per_module = pd.read_csv(OUTPUT_DIR / "counts_seeds_per_module.csv")
df_seeds_per_module.columns = ['KEGG_module', 'num_seeds']
plt.figure(figsize=(10,5))
sns.barplot(data=df_seeds_per_module.head(topN), x='KEGG_module', y='num_seeds', color='slateblue')
plt.xticks(rotation=45, ha='right')
plt.xlabel('KEGG Module')
plt.ylabel('Number of Seed Compounds')
plt.title(f'Top {topN} KEGG Modules by Number of Seed Compounds')
plt.tight_layout()
plt.show()

# Histogram: Genome coverage (% seeds/non-seeds per compound)
df_cov = pd.read_csv(OUTPUT_DIR / "genome_coverage_per_compound.csv")
plt.figure(figsize=(6,4))
df_cov['percent_as_seed'].hist(bins=20, alpha=0.6, label='Seeds')
df_cov['percent_as_non_seed'].hist(bins=20, alpha=0.6, label='Non-Seeds')
plt.xlabel('Percentage of Genomes')
plt.ylabel('Number of Compounds')
plt.legend()
plt.title('Distribution of Genome Coverage (%) for Seeds and Non-Seeds')
plt.tight_layout()
plt.show()

# Heatmap: Participation of top-N seeds in top-N KEGG pathways
matrix = pd.read_csv(OUTPUT_DIR / "seed_pathway_participation_matrix.csv", index_col=0)
topN_seeds = matrix.sum(axis=1).sort_values(ascending=False).head(topN).index
topN_pws = matrix.sum(axis=0).sort_values(ascending=False).head(topN).index
subset = matrix.loc[topN_seeds, topN_pws]
plt.figure(figsize=(12,8))
sns.heatmap(subset, cmap="Blues", cbar_kws={'label': 'Participation (1=present)'})
plt.xlabel('KEGG Pathway')
plt.ylabel('Seed Compound')
plt.title(f'Participation of Top {topN} Seeds in Top {topN} KEGG Pathways')
plt.tight_layout()
plt.show()
