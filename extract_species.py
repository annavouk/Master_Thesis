import pandas as pd
import matplotlib.pyplot as plt

#df = pd.read_csv('low_outliers_metadata.csv', low_memory=False)
df = pd.read_csv('high_outliers_metadata.csv', low_memory=False)

def extract_species(df):
    species = {
        part.split('__')[1]
        for taxonomy in df['gtdb_taxonomy'].dropna()
        for part in taxonomy.split(';')
        if part.startswith('s__')
    }

    return sorted(species)

species_list = extract_species(df)

species_df = pd.DataFrame(species_list, columns=['Species'])

print(species_df)

species_counts = species_df['Species'].value_counts()

# Plotting the bar chart
plt.figure(figsize=(12, 6))
species_counts.plot(kind='bar', color='skyblue')
plt.title('Species Distribution')
plt.xlabel('Species')
plt.ylabel('Frequency')
plt.xticks(rotation=90)  # Rotate the x-axis labels if they are long
plt.show()

# Plotting the pie chart
plt.figure(figsize=(8, 8))
species_counts.plot(kind='pie', autopct='%1.1f%%', startangle=90, cmap='Set3', legend=False)
plt.title('Species Distribution')
plt.ylabel('')  # Removing the y-axis label
plt.show()
