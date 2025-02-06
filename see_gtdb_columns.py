import pandas as pd

gtdb_df = pd.read_csv('gtdb_metadata_r207.tsv', sep='\t', low_memory=False)

column_names = gtdb_df.columns

print(column_names)
