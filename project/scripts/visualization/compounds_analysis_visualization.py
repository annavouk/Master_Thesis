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

vals = df_pathways_per_seed['num_kegg_pathways']
N = len(vals)

plt.figure(figsize=(6,4))
ax = plt.gca()
n, bins, patches = ax.hist(vals, bins=20, color='teal', edgecolor='black', log=True)
plt.xlabel('Number of KEGG Pathways per Seed')
plt.ylabel('Number of Seed Compounds (log scale)')
plt.yscale('log')
plt.title('Distribution of KEGG Pathways per Seed (N = {N} seed compounds)')
plt.tight_layout()

# Annotation
mean = vals.mean()
median = vals.median()
std = vals.std()
minv = vals.min()
maxv = vals.max()

plt.axvline(median, color='red', linestyle='--', label=f'Median = {median:.2f}')
plt.axvline(mean, color='purple', linestyle=':', label=f'Mean = {mean:.2f}')
plt.axvspan(mean-std, mean+std, color='purple', alpha=0.08, label=f'±1 STD = {std:.2f}')
plt.legend()
plt.annotate(f"Min: {minv:.0f}", xy=(minv, 0), xytext=(minv, 2), color='black', fontsize=10, rotation=90)
plt.annotate(f"Max: {maxv:.0f}", xy=(maxv, 0), xytext=(maxv, 2), color='black', fontsize=10, rotation=90)

plt.text(
    0.97, 0.95,
    f'N = {N}\nMean = {mean:.2f}\nMedian = {median:.2f}\nSTD = {std:.2f}',
    ha='right', va='top', transform=ax.transAxes,
    fontsize=12, bbox=dict(facecolor='white', edgecolor='gray', alpha=0.7)
)

# Add label (count) on top of each bar
for i in range(len(patches)):
    if n[i] > 0:
        plt.text(
            (bins[i]+bins[i+1])/2, n[i], f'{int(n[i])}',
            ha='left', va='bottom', fontsize=9, color='black'
        )

plt.tight_layout()
plt.show()

print("KEGG pathways per seed compound:")
print(f" N = {N}")
print(f" Mean = {mean:.2f}")
print(f" Median = {median:.2f}")
print(f" STD = {std:.2f}")
print(f" Min = {minv}")
print(f" Max = {maxv}")

# Barplot: Seeds per KEGG pathway (top 20)
df_seeds_per_pathway = pd.read_csv(OUTPUT_DIR / "counts_seeds_per_pathway.csv")
df_seeds_per_pathway.columns = ['KEGG_pathway', 'num_seeds']
topN = 20
N = len(df_seeds_per_pathway)

plt.figure(figsize=(10,5))
ax1 = sns.barplot(data=df_seeds_per_pathway.head(topN), x='KEGG_pathway', y='num_seeds', color='steelblue')
plt.xticks(rotation=45, ha='right')
plt.xlabel('KEGG Pathway')
plt.ylabel('Number of Seed Compounds')
plt.title(f'Top {topN} KEGG Pathways by Number of Seed Compounds (N = {N} pathways)')

# Label
for i, val in enumerate(df_seeds_per_pathway.head(topN)['num_seeds']):
    ax1.text(i, val + 0.5, str(val), ha='center', va='bottom', color='black', fontsize=9)

plt.tight_layout()
plt.show()

# Barplot: Seeds per KEGG module (top 20)
df_seeds_per_module = pd.read_csv(OUTPUT_DIR / "counts_seeds_per_module.csv")
df_seeds_per_module.columns = ['KEGG_module', 'num_seeds']

plt.figure(figsize=(10,5))
ax2 = sns.barplot(data=df_seeds_per_module.head(topN), x='KEGG_module', y='num_seeds', color='slateblue')
plt.xticks(rotation=45, ha='right')
plt.xlabel('KEGG Module')
plt.ylabel('Number of Seed Compounds')
plt.title(f'Top {topN} KEGG Modules by Number of Seed Compounds (N = {len(df_seeds_per_module)})')

# Label
for i, val in enumerate(df_seeds_per_module.head(topN)['num_seeds']):
    ax2.text(i, val + 0.5, str(val), ha='center', va='bottom', color='black', fontsize=9)

plt.tight_layout()
plt.show()

# Histogram: Genome coverage (% seeds/non-seeds per compound)
df_cov = pd.read_csv(OUTPUT_DIR / "genome_coverage_per_compound.csv")

plt.figure(figsize=(6,4))
ax = df_cov['percent_as_seed'].hist(bins=20, alpha=0.6, label='Seeds')
df_cov['percent_as_non_seed'].hist(bins=20, alpha=0.6, label='Non-Seeds')
plt.xlabel('Percentage of Genomes')
plt.ylabel('Number of Compounds')
plt.title(f'Distribution of Genome Coverage (%) for Seeds and Non-Seeds\n(Total compounds = {len(df_cov):,}, Total genomes = 33,755)')

mean_seed = df_cov['percent_as_seed'].mean()
mean_nonseed = df_cov['percent_as_non_seed'].mean()
plt.axvline(mean_seed, color='red', linestyle='--', label=f'Mean seed = {mean_seed:.1f}%')
plt.axvline(mean_nonseed, color='blue', linestyle='--', label=f'Mean non-seed = {mean_nonseed:.1f}%')
plt.legend()
plt.tight_layout()
plt.show()

print("Genome coverage (% of genomes):")
print("Seeds:")
print(df_cov['percent_as_seed'].describe())
print("Non-seeds:")
print(df_cov['percent_as_non_seed'].describe())

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
