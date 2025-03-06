import pickle
import csv
import requests
import json
import time
import pandas as pd
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

# Extract patric_ids list from pkl
pickle_file = 'seeds_binary_per_patric.pckl'
df = pd.read_pickle(pickle_file)

patric_ids_list = sorted(list(df.index))

# Metadata for each patric_id
output_file = "patric_ids_metadata.json"
failed_ids_file = "failed_patric_ids.txt"
metadata_url = "https://www.patricbrc.org/api/genome/"

patric_metadata = {}
failed_patric_ids = []

retry_strategy = Retry(
	total=5,
	backoff_factor=1,
	status_forcelist=[429, 500, 502, 503, 504],
	allowed_methods=["HEAD", "GET", "OPTIONS"]
)
adapter = HTTPAdapter(max_retries=retry_strategy)
http = requests.Session()
http.mount("https://", adapter)

def fetch_metadata(patric_id):
	try:
		response = http.get(f"{metadata_url}{patric_id}", timeout=10)
		if response.status_code == 200:
			return response.json()
		else:
			failed_patric_ids.append(patric_id)
			return None
	except requests.exceptions.RequestException as e:
		print(f"Error fetching data for {patric_id}: {e}")
		return None

for patric_id in patric_ids_list:
	metadata = fetch_metadata(patric_id)
	if metadata:
		patric_metadata[patric_id] = metadata
	time.sleep(1)

# Save output json and failed_patric_ids
with open(output_file, 'w') as jsonfile:
	json.dump(patric_metadata, jsonfile, indent=4)

if failed_patric_ids:
	with open(failed_ids_file, 'w') as f:
		for pid in failed_patric_ids:
			f.write(pid + "\n")
	print(f"Failed PATRIC IDs saved to {failed_ids_file}")

print(f"Metadata saved to {output_file}")
