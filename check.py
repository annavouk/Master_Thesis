import pandas as pd

df = pd.read_csv('gtdb_metadata_compact.csv')

df['values_equal'] = df['accession'] == df['ncbi_genbank_assembly_accession']

print(df[df['values_equal'] == False])
