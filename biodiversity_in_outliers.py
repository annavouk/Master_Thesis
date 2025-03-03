import pandas as pd
import matplotlib.pyplot as plt

#df = pd.read_csv("low_outliers_gtdb_taxonomy.csv", low_memory = False)
df = pd.read_csv("high_outliers_gtdb_taxonomy.csv", low_memory = False)

# Extract genera
def extract_genus(df):
	genus = [
		part.split('__')[1]
		for taxonomy in df['gtdb_taxonomy'].dropna()
		for part in taxonomy.split(';')
		if part.startswith('g__')
	]
	return genus


# Extract  classes
def extract_class(df):
	classes = [
		part.split('__')[1]
		for taxonomy in df['gtdb_taxonomy'].dropna()
		for part in taxonomy.split(';')
		if part.startswith('c__')
	]
	return classes

# Genus richness and most abundant genera
df['Genus'] = extract_genus(df)

genus_counts = df['Genus'].value_counts()

num_unique_genera = df['Genus'].nunique()
print(f"Total number of different Genera: {num_unique_genera}")

plt.figure(figsize=(12, 8))
genus_counts.head(10).plot(kind='bar', color='lightcoral')
plt.title('Top 10 Most Abundant Genera in high outliers')
plt.xlabel('Genus')
plt.ylabel('Count')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()

# Classes richness and most abundant classes
df['Classes'] = extract_class(df)

classes_counts = df['Classes'].value_counts()

num_unique_classes = df['Classes'].nunique()
print(f"Total number of different Classes: {num_unique_classes}")

plt.figure(figsize=(12, 8))
classes_counts.head(10).plot(kind='bar', color='skyblue')
plt.title('Top 10 Most Abundant Classes in high outliers')
plt.xlabel('Class')
plt.ylabel('Count')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
