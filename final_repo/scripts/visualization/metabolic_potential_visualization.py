"""
Visualize the metabolic potential of each genome.

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import ( ggplot, aes, element_text, labs, theme, geom_bar, geom_histogram, ggsave, geom_violin,
    theme, labs, element_text, scale_y_continuous)


# Load the data
df = pd.read_csv("../analysis/metabolic_potential_summary.csv", low_memory=False)


# Function to split and clean taxonomy
def split_and_clean_taxonomy(df, taxonomy_col_index):
    taxonomy_series = df.iloc[:, taxonomy_col_index].astype(str)
    taxonomy_split = taxonomy_series.str.split(';', expand=True)

    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxonomy_split.columns = levels[:taxonomy_split.shape[1]]

    taxonomy_cleaned = taxonomy_split.apply(lambda col: col.str.replace(r'^[a-z]__', '', regex=True))
    return taxonomy_cleaned

# Use the function
#taxonomy_df = split_and_clean_taxonomy(df, 5)

# Merge taxonomy back into main dataframe
#df = pd.concat([df, taxonomy_df], axis=1)

# Compute median seed/non-seed ratio per phylum
#phylum_ratio = df.groupby('phylum', as_index=False)['Ratio'].median()

# Sort for plotting
#phylum_ratio = phylum_ratio.sort_values(by='Ratio', ascending=False)

# Create barplot with plotnine
#p = (
#    ggplot(phylum_ratio, aes(x='reorder(phylum, -Ratio)', y='Ratio')) +
#    geom_bar(stat='identity', fill='salmon') +
#    labs(
#        title='Median Seed/Non-Seed Ratio Across Phyla',
#        x='Phylum',
#        y='Seed/Non-Seed Ratio'
#    ) +
#    theme(
#        axis_text_x=element_text(rotation=45, ha='right'),
#        figure_size=(20, 10)
#    )
#)

#p.save("median_ratio_per_phylum.png", dpi=300)

# Create violin plot with plotnine
#p1 = (
#    ggplot(df, aes(x='class', y='Ratio')) +
#    geom_violin(fill='lightblue') +
#    labs(
#        title='Distribution of Seed/Non-Seed Ratio Across Classes',
#        x='Class',
#        y='Seed/Non-Seed Ratio'
#    ) +
#    theme(
#        axis_text_x=element_text(rotation=45, ha='right', size=8),
#        figure_size=(14, 6)
#    ) +
#    scale_y_continuous(expand=(0.01, 0))  # Tighter y-axis
#)

# Save the figure
#p1.save("violin_ratio_per_classes.png", dpi=300)

#####################
### Visualization ###
#####################

# Histogram of Distribution of Total Seeds/nonSeeds/Ratio per Genome
plt.figure(figsize=(8, 5))
sns.histplot(df["Total_Seeds"], bins=30, kde=True)
plt.xlabel("Total Seeds per Genome")
plt.ylabel("Frequency")
plt.title("Distribution of Total Seeds per Genome")
plt.show()

# Boxplot Total Seeds/nonSeeds/Ratio per Genome
plt.figure(figsize=(6, 4))
sns.boxplot(x=df["Total_Seeds"])
plt.xlabel("Total Seeds per Genome")
plt.title("Boxplot Total Seeds per Genome")
plt.show()

#print(df["Total_Seeds"].describe())

