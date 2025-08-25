import pandas as pd

# Μικρό test DataFrame
data = {
    "cpd00002": [0, 1, 0],
    "cpd00007": [1, 0, 0],
    "cpd00015": [0, 0, 1],
}
index = ["100.11", "1000565.3", "1000566.3"]  # genome IDs

df = pd.DataFrame(data, index=index)
print(df)

# Save as pickle
df.to_pickle("test_patric.pkl")
