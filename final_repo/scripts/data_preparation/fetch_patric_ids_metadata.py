"""
Script to fetch metadata for PATRIC genome IDs using the PATRIC API.
Retrieves data from a pickle file, queries the API, and saves the results in JSON.
"""

import pickle
import requests
import json
import time
import pandas as pd
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

# Constants
PICKLE_FILE = "input/seeds_binary_per_patric.pckl" # Input
OUTPUT_FILE = "metadata/raw/patric_ids_metadata.json"
FAILED_IDS_FILE = "metadata/raw/failed_patric_ids.txt"
METADATA_URL = "https://www.patricbrc.org/api/genome/"

# Load PATRIC IDs
df = pd.read_pickle(PICKLE_FILE)
patric_ids_list = sorted(list(df.index))

# Containers for metadata and failures
patric_metadata = {}
failed_patric_ids = []

# Set up retry strategy for requests
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
    """
    Fetch metadata for a given PATRIC ID from the PATRIC API.

    Args:
        patric_id (str): The PATRIC genome ID.

    Returns:
        dict or None: Parsed JSON metadata if successful, None otherwise.
    """
    try:
        response = http.get(f"{METADATA_URL}{patric_id}", timeout=10)
        if response.status_code == 200:
            return response.json()
        else:
            failed_patric_ids.append(patric_id)
            return None
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data for {patric_id}: {e}")
        failed_patric_ids.append(patric_id)
        return None


# Fetch metadata for each ID
for patric_id in patric_ids_list:
    metadata = fetch_metadata(patric_id)
    if metadata:
        patric_metadata[patric_id] = metadata
    time.sleep(1)

# Save fetched metadata to JSON
with open(OUTPUT_FILE, 'w') as jsonfile:
    json.dump(patric_metadata, jsonfile, indent=4)

# Save failed IDs to text file
if failed_patric_ids:
    with open(FAILED_IDS_FILE, 'w') as f:
        for pid in failed_patric_ids:
            f.write(pid + "\n")
    print(f"Failed PATRIC IDs saved to {FAILED_IDS_FILE}")

print(f"Metadata saved to {OUTPUT_FILE}")

