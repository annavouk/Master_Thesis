import pandas as pd

df = pd.read_csv('gtdb_metadata_r207.tsv', sep='\t', low_memory=False)

accession_gtdb_taxonomy = df.iloc[:, [0, 16]]

accession_gtdb_taxonomy.columns = ['assembly_accession', 'gtdb_taxonomy']

accession_gtdb_taxonomy.loc[:, 'assembly_accession'] = accession_gtdb_taxonomy['assembly_accession'].str.replace('GB_|RS_', '', regex=True)

print(accession_gtdb_taxonomy.head())

accession_gtdb_taxonomy.to_csv('gtdb_accession_taxonomy.csv', index=False)
