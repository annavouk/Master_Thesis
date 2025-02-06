import json
import pandas as pd

with open('patric_ids_metadata.json', 'r') as file:
    data = json.load(file)

df = pd.DataFrame(data)

print(df.head())

