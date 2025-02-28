import pandas as pd
import numpy as np

patric_df = pd.read_csv('patric_to_assembly.csv')

output_df = pd.DataFrame(columns=['patric_id', 'refseq_accession', 'ncbi_genbank_accession'])

def fill_accession_type(row):
	refseq = np.nan
	genbank = np.nan

	if row['assembly_accession'].startswith('GCA_'):
		refseq = row['assembly_accession']
	elif row['assembly_accession'].startswith('GCF_'):
		genbank = row['assembly_accession']

	return pd.Series([row['patric_id'], refseq, genbank])

output_df = patric_df.apply(fill_accession_type, axis=1)

output_df.columns = ['patric_id', 'refseq_accession', 'ncbi_genbank_accession']

print(output_df)

output_df.to_csv('patric_refseq_genbank_accessions.csv', index=False)
