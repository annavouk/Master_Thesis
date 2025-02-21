import pandas as pd
import json

with open('patric_ids_metadata.json', 'r') as f:
    patric_data = json.load(f)

patric_df = pd.DataFrame.from_dict(patric_data, orient='index').reset_index()
patric_df.rename(columns={'index': 'patric_id'}, inplace=True)

gtdb_df = pd.read_csv('gtdb_accession_taxonomy.csv')

merged_df = pd.merge(patric_df, gtdb_df, on='assembly_accession', how='left')

print(merged_df.head())
merged_df.to_csv('merged_metadata.csv', index=False)
