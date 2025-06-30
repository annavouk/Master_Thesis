"""

"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from config import SEEDS_PICKLE, NON_SEEDS_PICKLE, COMPACT_METADATA, OUTPUT_DIR
from utils import load_data



df = pd.read_csv("/home/annavouk/master_thesis/project/output/genome_coverage_per_compound.csv")


top_n = 20
top_seeds = df.sort_values('percent_as_seed', ascending=False).head(top_n)

plt.figure(figsize=(10, 5))
sns.barplot(data=top_seeds, x='ModelSEED ID', y='percent_as_seed', color='teal')
plt.xlabel("Compound (ModelSEED ID)")
plt.ylabel("Percentage of Genomes as Seed")
plt.title(f"Top {top_n} Most Frequent Seed Compounds Across Genomes")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()

