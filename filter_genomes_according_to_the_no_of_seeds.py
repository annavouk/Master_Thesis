import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_pickle('seeds_to_non_seeds_df.pkl')
#print(df.head())

# Check distribution
# Histogram
plt.figure(figsize=(8, 5))
sns.histplot(df["Total_Seeds"], bins=30, kde=True)
plt.xlabel("Total Seeds per Genome")
plt.ylabel("Frequency")
plt.title("Distribution of Total Seeds per Genome")
plt.show()

# Boxplot
plt.figure(figsize=(6, 4))
sns.boxplot(x=df["Total_Seeds"])
plt.xlabel("Total Seeds per Genome")
plt.title("Boxplot Total Seeds per Genome")
plt.show()

print(df["Total_Seeds"].describe())

# Calculate low threshold (Q1) and high threshold (Q3)
df.sorted = df.sort_values(by="Total_Seeds")
Q1 = df.sorted["Total_Seeds"].quantile(0.25)
Q3 = df.sorted["Total_Seeds"].quantile(0.75)
print(f"Q1:{Q1}, Q3: {Q3}")

# Categorization of the genomes based on the no of seeds
# Category 1: genomes with no of seeds < Q1 or the least seeds (metabolic independence?/ genes involved in many metabolic pathways -> provide itself the needed compounds?)
# Category 2: genomes with no of seeds > Q3 or the most seeds (metabolic dependent?/ genes involved in a few metabolic pathways -> obtaining from the environment or from other microorganisms the needed compounds)
# Category 3: genomes with Q1 < no seeds < Q3 or moderate no of seeds
# Phylogenetic diversity in each category and between the categories?
df['category'] = df['Total_Seeds'].apply(lambda x: 1 if x < Q1 else (2 if x > Q3 else 3))

category_counts = df['category'].value_counts().sort_index()

print("Category distribution:", "/n", category_counts)

# Save each category to a file
category_1 = df[df['category'] == 1]
category_1.to_csv('category_1_genomes.csv', index=True)

category_2 = df[df['category'] == 2]
category_2.to_csv('category_2_genomes.csv', index=True)

category_3 = df[df['category'] == 3]
category_3.to_csv('category_3_genomes.csv', index=True)

# Check
print(category_1.head())
print(category_2.head())
print(category_3.head())
