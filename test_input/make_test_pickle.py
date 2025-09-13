import pandas as pd

# Seed test df
seed_data = {
    "cpd00002": [0, 1, 0],
    "cpd00007": [1, 0, 0],
    "cpd00015": [0, 0, 1],
}
index = ["100.11", "1000565.3", "1000566.3"]  # genome IDs

seed_df = pd.DataFrame(seed_data, index=index)
print(seed_df)

# Save as pickle
seed_df.to_pickle("test_seed_patric.pckl")

non_seed_data = {
    "cpd00002": [1, 0, 1],
    "cpd00007": [1, 1, 0],
    "cpd00008": [0, 0, 1],
}
index = ["100.11", "1000565.3", "1000566.3"]  # genome IDs

non_seed_df = pd.DataFrame(non_seed_data, index=index)
print(non_seed_df)

# Save as pickle
non_seed_df.to_pickle("test_non_seed_patric.pckl")
