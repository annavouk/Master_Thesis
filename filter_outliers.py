import pandas as pd

df1 = pd.read_csv('seeds_to_non_seeds.csv')
df1.rename(columns={'PATRIC': 'patric_id'}, inplace=True)
data = pd.read_csv("compact_metadata.csv", low_memory=False)

# Filter groups 0-5% and 95-100% of distribution
df1_sorted = df1.sort_values(by="Total_Seeds")
Q1 = df1_sorted["Total_Seeds"].quantile(0.05)
Q2 = df1_sorted["Total_Seeds"].quantile(0.95)

low_outliers = df1_sorted[df1_sorted["Total_Seeds"] <= Q1]
high_outliers = df1_sorted[df1_sorted["Total_Seeds"] >= Q2]

# Extract patric_id list
low_outliers_list = low_outliers["patric_id"].tolist()
high_outliers_list = high_outliers["patric_id"].tolist()

# Reset index correctly
low_outliers_reset = low_outliers.reset_index(drop=True)
high_outliers_reset = high_outliers.reset_index(drop=True)

# Ensure 'patric_id' is a string
def format_patric_id(x):
        try:
                return f"{float(x):.3f}"
        except ValueError:
                return str(x)

low_outliers_reset['patric_id'] = low_outliers_reset['patric_id'].astype(str).str.strip().apply(format_patric_id)
high_outliers_reset['patric_id'] = high_outliers_reset['patric_id'].astype(str).str.strip().apply(format_patric_id)
data['patric_id'] = data['patric_id'].astype(str).str.strip().apply(format_patric_id)

# Merge with metadata
low_outliers_gtdb_taxonomy = pd.merge(low_outliers_reset, data, on='patric_id', how='left').drop_duplicates(subset=['patric_id'])
high_outliers_gtdb_taxonomy = pd.merge(high_outliers_reset, data, on='patric_id', how='left').drop_duplicates(subset=['patric_id'])

# Save output
low_outliers_gtdb_taxonomy.to_csv('low_outliers_gtdb_taxonomy.csv', index=False)
high_outliers_gtdb_taxonomy.to_csv('high_outliers_gtdb_taxonomy.csv', index=False)

print(high_outliers_gtdb_taxonomy)
print(low_outliers_gtdb_taxonomy)
