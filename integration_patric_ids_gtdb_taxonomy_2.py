import pandas as pd

gtdb_df = pd.read_csv('gtdb_metadata_compact.csv')
patric_df = pd.read_csv('patric_to_assembly.csv')

# Remove the prefix (first 4 characters) from the accession columns in both dfs
patric_df['assembly_accession_numeric'] = patric_df['assembly_accession'].str.slice(4)
gtdb_df['refseq_accession_numeric'] = gtdb_df['refseq_accession'].str.slice(4)

merged_df = pd.merge(patric_df, gtdb_df, left_on='assembly_accession_numeric', right_on='refseq_accession_numeric', how='left')

merged_df.to_csv('genome_size_merged_metadata.csv', index=False)

print(merged_df)
