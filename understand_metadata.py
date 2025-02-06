import pandas as pd

metadata = pd.read_csv('merged_metadata.csv', low_memory=False)

print(metadata.head())

print(metadata.columns)
