import pandas as pd
import pickle
import json
from pathlib import Path


def load_data(filepath, filetype=None, dtype=None):
    """
    Load data from various formats: csv, tsv, pickle, json.
    If filetype is None, infer from file extension.
    """
    filepath = Path(filepath)

    if not filetype:
        suffix = filepath.suffix if isinstance(filepath, Path) else Path(filepath).suffix
        if suffix == '.csv':
            filetype = 'csv'
        elif suffix == '.tsv' or suffix == '.pathways':
            filetype = 'tsv'
        elif suffix == '.pckl'or suffix == '.pkl':
            filetype = 'pickle'
        elif suffix == '.json':
            filetype = 'json'
        else:
            raise ValueError(f"Cannot infer filetype from extension for file: {filepath}")

    if filetype == 'csv':
        return pd.read_csv(filepath, low_memory=False, dtype=dtype)
    elif filetype == 'tsv':
        return pd.read_csv(filepath, sep='\t', low_memory=False, dtype=dtype)
    elif filetype == 'pickle':
        with open(filepath, 'rb') as f:
            return pickle.load(f)
    elif filetype == 'json':
        with open(filepath, 'r') as f:
            return json.load(f)
    else:
        raise ValueError(f"Unsupported filetype: {filetype}")
