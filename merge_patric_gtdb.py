import pandas as pd

patric = pd.read_csv("patric_refseq_genbank_accessions.csv", low_memory=False)

gtdb = pd.read_csv("gtdb_metadata_compact.csv", low_memory=False)

merged_data = pd.merge(patric, gtdb, on='refseq_accession', how='inner')
merged_data = pd.merge(patric, gtdb, on='ncbi_genbank_accession', how='inner')

print(merged_data)

merged_data.to_csv("patric_ids_with_gtdb_taxonomy_correct.csv", index=False)
