import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('low_outliers_gtdb_taxonomy.csv', low_memory=False)

df[['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = df['gtdb_taxonomy'].str.split(';', expand=True)

# Scatter plot: genome_size vs. Total_Seeds
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x='genome_size_x', y='Total_Seeds', hue='Genus', size='Genus', sizes=(20, 200))
plt.title(' Genome size vs. Total number of Seeds')
plt.xlabel('Genome Size')
plt.ylabel('Total number of Seeds')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()

# Scatter plot: genome_size vs. Total_non_Seeds
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x='genome_size_x', y='Total_non_Seeds', hue='Genus', size='Genus', sizes=(20, 200))
plt.title(' Genome size vs. Total number of non Seeds')
plt.xlabel('Genome Size')
plt.ylabel('Total number of non Seeds')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
