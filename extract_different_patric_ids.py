
import pandas as pd

df = pd.read_csv("missing_taxonomy.csv", low_memory=False)
sorted_df = df.sort_values(by='patric_id')
patric_ids = sorted_df["patric_id"]

# Save to a new CSV file
patric_ids.to_csv("different_patric_ids.csv", index=False)

# Print the first few IDs
print(patric_ids.head())
