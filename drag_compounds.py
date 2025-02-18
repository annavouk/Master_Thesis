import pickle
import pandas as pd

pickle_file_1 = 'seeds_binary_per_patric.pckl'
df_1 = pd.read_pickle(pickle_file_1)

pickle_file_2 = 'non_seeds_binary_per_patric.pckl'
df_2 = pd.read_pickle(pickle_file_2)

unique_seed_set_list = list(set(df_1.columns.tolist() + df_2.columns.tolist()))

unique_df = pd.DataFrame(unique_seed_set_list, columns=['cpd_id'])

unique_df.to_csv('unique_compounds.csv', index=False, header=True)

print(unique_df.head())
print(unique_df.tail())
