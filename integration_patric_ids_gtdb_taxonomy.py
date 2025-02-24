import pandas as pd
import json

with open('patric_ids_metadata.json', 'r') as f:
    patric_data = json.load(f)

patric_df = pd.DataFrame.from_dict(patric_data, orient='index').reset_index()
patric_df.rename(columns={'index': 'patric_id'}, inplace=True)

gtdb_df = pd.read_csv('patric_ids_with_gtdb_taxonomy.csv')

patric_df['patric_id'] = patric_df['patric_id'].astype(str)
gtdb_df['patric_id'] = gtdb_df['patric_id'].astype(str)

merged_df = pd.merge(patric_df, gtdb_df, on='patric_id', how='left')

print(merged_df.head())

merged_df.to_csv('merged_metadata.csv', index=False)
