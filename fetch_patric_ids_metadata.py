import csv
import requests
import json
import time
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

patric_ids_file = "patric_ids.csv"
output_file = "patric_ids_metadata.json"
metadata_url = "https://www.patricbrc.org/api/genome/"

patric_metadata = {}

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
            print(f"Failed to retrieve metadata for {patric_id}, status code: {response.status_code}")
            return None
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data for {patric_id}: {e}")
        return None

with open(patric_ids_file, 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        patric_id = row[0]
        metadata = fetch_metadata(patric_id)
        if metadata:
            patric_metadata[patric_id] = metadata
        time.sleep(1)

with open(output_file, 'w') as jsonfile:
    json.dump(patric_metadata, jsonfile, indent=4)

print(f"Metadata saved to {output_file}")
