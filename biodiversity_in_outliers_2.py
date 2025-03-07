import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('seeds_to_non_seeds_genome_size.csv', low_memory=False)

df[['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = df['gtdb_taxonomy'].str.split(';', expand=True)

# Scatter plot: genome_size vs. Total_Seeds in low outliers
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x='genome_size', y='Total_Seeds', hue='Genus', size='Class', sizes=(20, 200))
plt.title(' Genome size vs. Total number of Seeds')
plt.xlabel('Genome Size')
plt.ylabel('Total number of Seeds')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

# Scatter plot: genome_size vs. Total_non_Seeds
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x='genome_size', y='Total_non_Seeds', hue='Genus', size='Class', sizes=(20, 200))
plt.title(' Genome size vs. Total number of non Seeds')
plt.xlabel('Genome Size')
plt.ylabel('Total number of non Seeds')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
