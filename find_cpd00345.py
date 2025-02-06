import pickle
import pandas as pd

pickle_file = 'seeds_binary_per_patric.pckl'

with open(pickle_file, 'rb') as f:
        data = pickle.load(f)

if 'cpd00345' in data.columns:
	print(data['cpd00345'])
else:
	print("Column 'cpd00345' not found.")
