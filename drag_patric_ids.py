import pickle
import pandas as pd

pickle_file = 'seeds_binary_per_patric.pckl'
df = pd.read_pickle(pickle_file)

patric_ids_list = list(df.index)

patric_ids_list = sorted(df.index)

pd.Series(patric_ids_list).to_csv('patric_ids.csv', index=False)
