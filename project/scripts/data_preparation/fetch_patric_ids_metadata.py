"""
Script to fetch metadata for PATRIC genome IDs using the PATRIC API.
Retrieves data from a pickle file, queries the API, and saves the results in JSON.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pickle
import requests
import json
import time
import pandas as pd
from config import SEEDS_PICKLE, PATRIC_METADATA_JSON, FAILED_PATRIC_IDS_TXT
from utils import load_data
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


# Constants
METADATA_URL = "https://www.patricbrc.org/api/genome/"

def setup_requests_session():
    """Set up a requests session with retry strategy"""
    retry_strategy = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["HEAD", "GET", "OPTIONS"]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session()
    http.mount("https://", adapter)
    return http

def fetch_metadata(patric_id, http):
    """Fetch metadata for a single PATRIC ID, with error handling."""
    try:
        response = http.get(f"{METADATA_URL}{patric_id}", timeout=10)
        if response.status_code == 200:
            return response.json()
        else:
            return None
    except requests.exceptions.RequestException as e:
        print(f"Request error for {patric_id}: {e}")
        return None

def fetch_all_metadata(patric_ids, http, delay=1):
    """Loop through all IDs, fetch metadata, and track failures."""
    patric_metadata = {}
    failed_patric_ids = []
    for i, patric_id in enumerate(patric_ids, 1):
        metadata = fetch_metadata(patric_id, http)
        if metadata:
            patric_metadata[patric_id] = metadata
        else:
            failed_patric_ids.append(patric_id)
        print(f"[{i}/{len(patric_ids)}] Done: {patric_id}", end="\r")
        time.sleep(delay)
    print()
    return patric_metadata, failed_patric_ids

def save_json(data, outpath):
    with outpath.open('w') as f:
        json.dump(data, f, indent=4)

def save_failed_ids(failed_ids, outpath):
    if failed_ids:
        with outpath.open('w') as f:
            for pid in failed_ids:
                f.write(pid + "\n")
        print(f"Failed PATRIC IDs saved to {outpath}")
    else:
        print("No failed PATRIC IDs.")

def main():
    http = setup_requests_session()
    df = load_data(SEEDS_PICKLE)
    patric_ids = sorted(list(df.index))
    patric_metadata, failed_patric_ids = fetch_all_metadata(patric_ids, http)
    save_json(patric_metadata, PATRIC_METADATA_JSON)
    save_failed_ids(failed_patric_ids, FAILED_PATRIC_IDS_TXT)
    print(f"Metadata saved to {PATRIC_METADATA_JSON}")

if __name__ == "__main__":
    main()