import pandas as pd

df = pd.read_csv('merged_metadata.csv', low_memory=False)

print(df.head())

middle_index = len(df) // 2
print(df.iloc[middle_index - 2 : middle_index + 3])

print(df.tail())
