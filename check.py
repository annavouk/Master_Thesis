import pandas as pd

df = pd.read_csv('gtdb_metadata_compact.csv')

# Remove the 'GCA' or 'GCF' prefixes and compare only the numeric part of the accession IDs
df['refseq_accession_base'] = df['refseq_accession'].str.replace(r'^GCA_|^GCF_', '', regex=True)
df['ncbi_genbank_accession_base'] = df['ncbi_genbank_accession'].str.replace(r'^GCA_|^GCF_', '', regex=True)

df['values_equal'] = df['refseq_accession_base'] == df['ncbi_genbank_accession_base']

print(df[df['values_equal'] == False])
