


## Building the binary seeds `.pickle` files

Starting from the `updated_seedsets_of_interest.pckl` and the `updated_non_seedsets_of_interest.pckl` files you may get from [Zenodo](https://doi.org/10.5281/zenodo.10562677), we had:

```python
import pickle
import pandas as pd

with open("updated_non_seedsets_of_interest.pckl", "rb") as f:
    df = pickle.load(f)

df = df.reset_index()
df.columns = ['PATRIC', 'NonSeedSet']

df_exploded = df.explode('NonSeedSet')

# pd.crosstab() computes a frequency table of the factors
# https://pandas.pydata.org/docs/reference/api/pandas.crosstab.html
df_binary = pd.crosstab(df_exploded['PATRIC'], df_exploded['NonSeedSet'])

# Export 
with open("non_seeds_binary_per_patric.pckl", "wb") as f:
    pickle.dump(df_binary ,f)

```


