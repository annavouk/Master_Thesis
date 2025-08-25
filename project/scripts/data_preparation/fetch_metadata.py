"""
Generic fetcher for metadata from KEGG, PATRIC, and PREGO APIs.

- Supports multiple data sources (KEGG compounds, PATRIC genomes, PREGO environments)
- Implements retry/caching and error logging
- Saves JSON outputs and lists of failed IDs

Usage:
    python3 fetch_metadata.py --source kegg
    python3 fetch_metadata.py --source patric
    python3 fetch_metadata.py --source prego
    python3 fetch_metadata.py --source kegg --out_json results/kegg.json --out_failed results/kegg_failed.txt
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import time
import json
import argparse
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

from config import (
    METADATA_DIR,
    CPDs_KEGG_DATASET_TSV,
    SEEDS_PICKLE,
    METABOLIC_POTENTIAL_1,
    PATRIC_METADATA_JSON,
    FAILED_PATRIC_IDS_TXT,
    PREGO_ENVIRONMENTS_JSON,
    FAILED_PREGO_IDS_TXT,
)
from utils import load_data

# Default outputs
KEGG_DATA_JSON = METADATA_DIR / "kegg_flat_entries.json"
FAILED_KEGG_IDS_TXT = METADATA_DIR / "failed_kegg_ids.txt"


# ------------------------
# Session + helpers
# ------------------------
def setup_requests_session():
    """Set up a requests session with retry strategy."""
    retry_strategy = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["HEAD", "GET", "OPTIONS"],
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session()
    http.mount("https://", adapter)
    return http


def save_json(data, outpath):
    """Save dictionary to JSON file with string keys."""
    data_str_keys = {str(k): v for k, v in data.items()}
    with open(outpath, "w") as f:
        json.dump(data_str_keys, f, indent=2)


def save_failed_ids(failed, outpath):
    """Save list of failed IDs to text file."""
    if failed:
        with open(outpath, "w") as f:
            for fid in failed:
                f.write(str(fid) + "\n")
        print(f"Failed IDs saved to {outpath}")
    else:
        print("No failed IDs.")


# ------------------------
# Parsers
# ------------------------
def parse_kegg_flat_text(text: str) -> dict:
    """Parse flat KEGG compound entry text to structured dict."""
    data = {}
    lines = text.split("\n")
    current_key = None
    for line in lines:
        if not line.strip():
            continue
        if line[:12].strip():
            current_key = line[:12].strip()
            value = line[12:].strip()
            data.setdefault(current_key, []).append(value)
        elif current_key:
            data[current_key].append(line[12:].strip())
    return data


def parse_prego_response(taxid: int, http, cache_dir: Path):
    """Fetch PREGO envs for one TaxID (literature + metagenomic)."""
    out = {"lit_envs": [], "meta_envs": []}
    endpoints = {
        "lit_envs": f"https://prego.hcmr.gr/Textmining?type1=-2&type2=-27&limit=100000&format=json&id1={taxid}",
        "meta_envs": f"https://prego.hcmr.gr/Experiments?type1=-2&type2=-27&limit=100000&format=json&id1={taxid}",
    }
    for kind, url in endpoints.items():
        cache_path = cache_dir / f"{taxid}_{kind}.json"
        if cache_path.exists():
            data = json.loads(cache_path.read_text())
        else:
            try:
                r = http.get(url, timeout=30)
                r.raise_for_status()
                cache_path.write_text(r.text)
                data = r.json()
            except Exception as e:
                print(f"Error fetching {url}: {e}")
                continue

        envs = set()
        if isinstance(data, list):
            records = [rec for rec in data if isinstance(rec, dict)]
        elif isinstance(data, dict) and "results" in data:
            records = data["results"]
        else:
            records = []

        for rec in records:
            if isinstance(rec, dict):
                for _, env_info in rec.items():
                    label = env_info.get("name", "").strip()
                    if label:
                        envs.add(label.lower())
        out[kind] = sorted(envs)
        time.sleep(0.01)
    return out


# ------------------------
# Generic fetch loop
# ------------------------
def fetch_one(http, url_template, id_, parser=None):
    """Fetch metadata for a single ID from a source."""
    try:
        response = http.get(url_template.format(id=id_), timeout=10)
        if response.status_code == 200:
            return (
                parser(response.text if parser == parse_kegg_flat_text else response)
                if parser
                else response.json()
            )
        else:
            print(f"Failed to fetch {id_}: HTTP {response.status_code}")
            return None
    except requests.exceptions.RequestException as e:
        print(f"Request error for {id_}: {e}")
        return None


def fetch_all(http, ids, url_template, parser=None, delay=0.01):
    """Fetch metadata for multiple IDs from a source."""
    data, failed = {}, []
    for i, id_ in enumerate(ids, 1):
        result = fetch_one(http, url_template, id_, parser=parser)
        if result:
            data[id_] = result
        else:
            failed.append(id_)
        print(f"[{i}/{len(ids)}] Done: {id_}", end="\r")
        time.sleep(delay)
    print()
    return data, failed


# ------------------------
# Sources config
# ------------------------
def load_kegg_ids():
    """Load KEGG compound IDs from the dataset TSV."""
    df = load_data(CPDs_KEGG_DATASET_TSV, filetype="tsv")
    df = df[df["KEGG_ID"].notna()]
    if df["KEGG_ID"].apply(lambda x: isinstance(x, str) and "," in x).any():
        df["KEGG_ID"] = df["KEGG_ID"].apply(lambda x: [k.strip() for k in x.split(",")])
    if df["KEGG_ID"].apply(lambda x: isinstance(x, list)).any():
        df = df.explode("KEGG_ID")
    df["KEGG_ID"] = df["KEGG_ID"].str.strip()
    df = df[df["KEGG_ID"].notna() & (df["KEGG_ID"] != "")]
    return df["KEGG_ID"].drop_duplicates().tolist()


def load_patric_ids():
    """Load PATRIC IDs from the seeds dataset."""
    df = load_data(SEEDS_PICKLE)
    return sorted(list(df.index))


def load_taxids():
    """Load unique NCBI taxids from the metabolic potential dataset."""
    df = load_data(METABOLIC_POTENTIAL_1, filetype="csv")
    return sorted(df["ncbi_taxid"].dropna().unique())


SOURCES = {
    "kegg": {
        "url_template": "https://rest.kegg.jp/get/compound:{id}",
        "parser": parse_kegg_flat_text,
        "id_loader": load_kegg_ids,
        "out_json": KEGG_DATA_JSON,
        "out_failed": FAILED_KEGG_IDS_TXT,
    },
    "patric": {
        "url_template": "https://www.patricbrc.org/api/genome/{id}",
        "parser": lambda r: r.json(),
        "id_loader": load_patric_ids,
        "out_json": PATRIC_METADATA_JSON,
        "out_failed": FAILED_PATRIC_IDS_TXT,
    },
    "prego": {
        "url_template": None,  # special handling
        "parser": None,
        "id_loader": load_taxids,
        "out_json": PREGO_ENVIRONMENTS_JSON,
        "out_failed": FAILED_PREGO_IDS_TXT,
    },
}


# ------------------------
# Main
# ------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Fetch metadata from KEGG, PATRIC or PREGO"
    )
    parser.add_argument(
        "--source",
        choices=["kegg", "patric", "prego"],
        required=True,
        help="Source database to fetch from",
    )
    parser.add_argument("--out_json", type=Path, help="Custom output JSON path")
    parser.add_argument("--out_failed", type=Path, help="Custom failed IDs path")
    args = parser.parse_args()

    cfg = SOURCES[args.source]
    out_json = args.out_json or cfg["out_json"]
    out_failed = args.out_failed or cfg["out_failed"]

    http = setup_requests_session()
    ids = cfg["id_loader"]()
    print(f"Found {len(ids)} IDs for {args.source}")

    if args.source == "prego":
        cache_dir = METADATA_DIR / "prego_cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        results, failed = {}, []

        for i, tx in enumerate(ids, 1):
            lit_file = cache_dir / f"{tx}_lit_envs.json"
            meta_file = cache_dir / f"{tx}_meta_envs.json"

            # If both cache files already exist -> read directly (no API request)
            if lit_file.exists() and meta_file.exists():
                envs = parse_prego_response(tx, http, cache_dir)
                results[tx] = envs
            else:
                # Otherwise fetch from API and store
                envs = parse_prego_response(tx, http, cache_dir)
                if not envs["lit_envs"] and not envs["meta_envs"]:
                    failed.append(tx)
                results[tx] = envs
                time.sleep(0.05)  # small delay only for real API requests

            print(f"[{i}/{len(ids)}] Done: {tx}", end="\r")

        print()

    else:
        results, failed = fetch_all(
            http, ids, cfg["url_template"], parser=cfg["parser"]
        )

    save_json(results, out_json)
    save_failed_ids(failed, out_failed)
    print(f"{args.source.upper()} data saved to {out_json}")


if __name__ == "__main__":
    main()
