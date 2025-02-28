import pandas as pd

patric_df = pd.read_csv('patric_to_assembly.csv')

output_df = pd.DataFrame(columns=['patric_id', 'refseq_accession', 'ncbi_genbank_accession'])

def fill_accession_type(row):
	refseq = row['assembly_accession']
	genbank = row['assembly_accession']
	return pd.Series([row['patric_id'], refseq, genbank])

output_df = patric_df.apply(fill_accession_type, axis=1)

output_df.columns = ['patric_id', 'refseq_accession', 'ncbi_genbank_accession']

print(output_df)

output_df.to_csv('patric_refseq_genbank_accessions.csv', index=False)
