import pandas as pd

patric_df = pd.read_csv('gtdb2patricIds.tsv', sep='\t', header=None)

patric_df.columns = ['assembly_accession', 'patric_id']

patric_df = patric_df.reset_index(drop=True)

gtdb_df = pd.read_csv('gtdb_accession_taxonomy.csv')

merged_df = pd.merge(patric_df, gtdb_df, on='assembly_accession', how='left')

print(merged_df)
merged_df.to_csv('patric_ids_with_gtdb_taxonomy.csv', index=False)

