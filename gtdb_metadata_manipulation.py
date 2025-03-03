import pandas as pd

df = pd.read_csv('gtdb_metadata_r207.tsv', sep='\t', low_memory=False)

accession_gtdb_taxonomy = df.iloc[:, [0, 13, 16, 54]]

accession_gtdb_taxonomy.columns = ['accession', 'genome_size', 'gtdb_taxonomy', 'ncbi_genbank_assembly_accession']

accession_gtdb_taxonomy.loc[:, 'accession'] = accession_gtdb_taxonomy['accession'].str.replace('GB_|RS_', '', regex=True)

accession_gtdb_taxonomy = accession_gtdb_taxonomy.rename(columns={'accession': "refseq_accession", 'genome_size': "genome_size", 'gtdb_taxonomy': "gtdb_taxonomy", 'ncbi_genbank_assembly_accession': "ncbi_genbank_accession"})

print(accession_gtdb_taxonomy.head())

accession_gtdb_taxonomy.to_csv('gtdb_metadata_compact.csv', index=False)
