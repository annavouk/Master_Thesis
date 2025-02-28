import pandas as pd

df1 = pd.read_pickle('seeds_to_non_seeds_df.pkl')
data = pd.read_csv("patric_ids_with_gtdb_taxonomy_correct.csv", low_memory=False) # patric_id is the first column

# Filter groups 0-5% and 95-100% of distribution
df1_sorted = df1.sort_values(by="Total_Seeds")
Q1 = df1_sorted["Total_Seeds"].quantile(0.05)
Q2 = df1_sorted["Total_Seeds"].quantile(0.95)

low_outliers = df1_sorted[df1_sorted["Total_Seeds"] <= Q1] # PATRIC is index
high_outliers = df1_sorted[df1_sorted["Total_Seeds"] >= Q2] # PATRIC is index

# To check later if they are all present
low_outliers_list = low_outliers.index.tolist()
high_outliers_list = high_outliers.index.tolist()

# Reset index as the first column and rename it to patric_id
low_outliers_reset = low_outliers.reset_index()
low_outliers_reset.rename(columns={low_outliers_reset.columns[0]: 'patric_id'}, inplace=True)

high_outliers_reset = high_outliers.reset_index()
high_outliers_reset.rename(columns={high_outliers_reset.columns[0]: 'patric_id'}, inplace=True)

# Convert 'patric_id' to string (object type) in both DataFrames to overcome missing 0s in IDs
low_outliers_reset['patric_id'] = low_outliers_reset['patric_id'].astype(str).str.strip()
high_outliers_reset['patric_id'] = high_outliers_reset['patric_id'].astype(str).str.strip()

data['patric_id'] = data['patric_id'].astype(str).str.strip()

low_outliers_gtdb_taxonomy = pd.merge(low_outliers_reset, data, on='patric_id', how='left')

low_outliers_gtdb_taxonomy = low_outliers_gtdb_taxonomy.drop_duplicates(subset=['patric_id'], keep='first')

low_outliers_gtdb_taxonomy.to_csv('low_outliers_gtdb_taxonomy.csv', index=False)

high_outliers_gtdb_taxonomy = pd.merge(high_outliers_reset, data, on='patric_id', how='left')

high_outliers_gtdb_taxonomy = high_outliers_gtdb_taxonomy.drop_duplicates(subset=['patric_id'], keep='first')

high_outliers_gtdb_taxonomy.to_csv('high_outliers_gtdb_taxonomy.csv', index=False)

# Check if all patric_ids of interest are present in the merged file
present_ids_low = low_outliers_gtdb_taxonomy['patric_id'].isin(low_outliers_list)

for patric_id, is_present in zip(low_outliers_list, present_ids_low):
	if not is_present:
		print(f"{patric_id} is NOT present in the merged file.")

present_ids_high = high_outliers_gtdb_taxonomy['patric_id'].isin(high_outliers_list)

for patric_id, is_present in zip(high_outliers_list, present_ids_high):
        if not is_present:
                print(f"{patric_id} is NOT present in the merged file.")

print(high_outliers_gtdb_taxonomy)
print(low_outliers_gtdb_taxonomy)
