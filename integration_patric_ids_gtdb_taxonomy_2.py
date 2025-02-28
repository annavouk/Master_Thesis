import pandas as pd
import json

gtdb_df = pd.read_csv('gtdb_metadata_compact.csv')

patric_df['patric_id'] = patric_df['patric_id'].astype(str)
gtdb_df['patric_id'] = gtdb_df['patric_id'].astype(str)

merged_df = pd.merge(patric_df, gtdb_df, on='patric_id', how='left')

print(merged_df)

merged_df.to_csv('merged_metadata.csv', index=False)
