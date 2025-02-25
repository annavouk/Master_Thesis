import pandas as pd

df = pd.read_csv('gtdb_metadata_r207.tsv', sep='\t', low_memory=False)

accession_gtdb_taxonomy = df.iloc[:, [0, 16, 54]]

accession_gtdb_taxonomy.columns = ['accession', 'gtdb_taxonomy', 'ncbi_genbank_assembly_accession']

accession_gtdb_taxonomy.loc[:, 'accession'] = accession_gtdb_taxonomy['accession'].str.replace('GB_|RS_', '', regex=True)

print(accession_gtdb_taxonomy.head())

accession_gtdb_taxonomy.to_csv('gtdb_metadata_compact.csv', index=False)
