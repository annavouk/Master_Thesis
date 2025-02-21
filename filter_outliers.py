import pandas as pd

df1 = pd.read_pickle('seeds_to_non_seeds_df.pkl')
metadata = pd.read_csv("merged_metadata.csv", low_memory=False) # patric_id is the first column
metadata = metadata.dropna(subset=['patric_id'])

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

# Convert 'patric_id' to string (object type) in both DataFrames
low_outliers_reset['patric_id'] = low_outliers_reset['patric_id'].astype(str)
metadata['patric_id'] = metadata['patric_id'].astype(str)

low_outliers_metadata = pd.merge(low_outliers_reset, metadata, on='patric_id', how='left')

low_outliers_metadata.to_csv('low_outliers_metadata.csv', index=False)

high_outliers_reset['patric_id'] = high_outliers_reset['patric_id'].astype(str)
metadata['patric_id'] = metadata['patric_id'].astype(str)

high_outliers_metadata = pd.merge(high_outliers_reset, metadata, on='patric_id', how='left')

high_outliers_metadata.to_csv('high_outliers_metadata.csv', index=False)

# Check if all patric_ids of interest are present in the merged file
present_ids_low = low_outliers_metadata['patric_id'].isin(low_outliers_list)

for patric_id, is_present in zip(low_outliers_list, present_ids_low):
	if not is_present:
		print(f"{patric_id} is NOT present in the metadata file.")

present_ids_high = high_outliers_metadata['patric_id'].isin(high_outliers_list)

for patric_id, is_present in zip(low_outliers_list, present_ids_high):
        if not is_present:
                print(f"{patric_id} is NOT present in the metadata file.")

