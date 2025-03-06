import pandas as pd
import json

with open('patric_ids_metadata.json', 'r') as f:
	patric_data = json.load(f)

patric_df = pd.DataFrame.from_dict(patric_data, orient='index').reset_index()
patric_df.rename(columns={'index': 'patric_id'}, inplace=True)
patric_df_filtered = patric_df[['patric_id', 'assembly_accession', 'genome_length']]
patric_df_filtered.loc[:, 'assembly_accession'] = patric_df_filtered['assembly_accession'].str.replace(r'^(GCA_|GCF_)', '', regex=True)

gtdb_df = pd.read_csv('gtdb_metadata_r207.tsv', sep='\t', low_memory=False)

accession_gtdb_taxonomy = gtdb_df[['accession', 'genome_size', 'gtdb_taxonomy', 'ncbi_genbank_assembly_accession']].copy()
accession_gtdb_taxonomy['accession'] = accession_gtdb_taxonomy['accession'].str.replace('GB_|RS_', '', regex=True)
accession_gtdb_taxonomy.rename(columns={'accession': 'assembly_accession'}, inplace=True)
accession_gtdb_taxonomy.loc[:, 'assembly_accession'] = accession_gtdb_taxonomy['assembly_accession'].str.replace(r'^(GCA_|GCF_)', '', regex=True)

merged_df = pd.merge(patric_df_filtered, accession_gtdb_taxonomy, on='assembly_accession', how='left')

merged_df.to_csv('compact_metadata.csv', index=False)
