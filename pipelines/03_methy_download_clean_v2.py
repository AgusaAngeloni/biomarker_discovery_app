#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_methy_download_clean_v2.py

Description
-----------
Downloads and filters TCGA HumanMethylation450K beta-value data.

The methylation matrix is filtered using:

    - CpG probes present in the cleaned manifest
    - Samples present in the cleaned phenotype table
    - CpG sites with at least 20% valid beta-values

Therefore, CpGs with more than 80% missing values are removed.

Processing is performed in chunks to minimize memory consumption.

Output
------
data/raw/methy.parquet
"""

# ============================================================
# Config
# ============================================================

from pathlib import Path
import requests
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


# ============================================================
# Data Dir
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = ROOT / "data" / "raw"
DATA_DIR.mkdir(parents=True, exist_ok=True)

METHYLATION_URL = (
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/"
    "GDC-PANCAN.methylation450.tsv"
)

TSV_PATH = DATA_DIR / "methylation_450k.tsv"
PARQUET_PATH = DATA_DIR / "methy.parquet"

MANIFEST_PATH = DATA_DIR / "manifest_clean.parquet"
PHENO_PATH = DATA_DIR / "pheno_clean.parquet"

CHUNK_SIZE = 5000

# Keeps CpGs with at least 20% valid data.
# Removes CpGs with more than 80% missing values.
MIN_VALID_FRACTION = 0.20

# ============================================================
# Functions
# ============================================================

def download_file():
    """
    Download methylation TSV file if it does not already exist.
    """
    if TSV_PATH.exists():
        print("✔ File already downloaded")
        return

    print("Downloading .tsv file")

    r = requests.get(METHYLATION_URL, stream=True)
    r.raise_for_status()

    with open(TSV_PATH, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)


def load_filters():
    """
    Load valid CpG sites and valid sample IDs.
    """
    manifest = pd.read_parquet(MANIFEST_PATH)
    pheno = pd.read_parquet(PHENO_PATH)

    valid_sites = set(manifest["site_id"])
    valid_samples = set(pheno["sample_id"])

    return valid_sites, valid_samples


def filter_missing_values(chunk):
    """
    Remove CpGs with more than 80% missing beta-values.

    Parameters
    ----------
    chunk : pd.DataFrame
        Methylation chunk with site_id and sample columns.

    Returns
    -------
    pd.DataFrame
        Filtered chunk.
    """

    sample_cols = [
        c for c in chunk.columns
        if c != "site_id"
    ]

    if len(sample_cols) == 0:
        raise ValueError(
            "No sample columns found after filtering by phenotype samples."
        )

    # Convert beta-values to numeric.
    # Non-numeric values are converted to NaN.
    chunk[sample_cols] = (
        chunk[sample_cols]
        .apply(pd.to_numeric, errors="coerce")
        .astype("float32")
    )

    n_samples = len(sample_cols)

    min_valid = max(
        1,
        int(np.ceil(n_samples * MIN_VALID_FRACTION))
    )

    n_valid = chunk[sample_cols].notna().sum(axis=1)

    chunk = chunk[
        n_valid >= min_valid
    ].copy()

    return chunk


def process_methylation(valid_sites, valid_samples):
    """
    Process methylation file by chunks and save filtered parquet.
    """

    writer = None

    total_input = 0
    total_after_manifest = 0
    total_after_missing_filter = 0

    for i, chunk in enumerate(
        pd.read_csv(
            TSV_PATH,
            sep="\t",
            chunksize=CHUNK_SIZE
        )
    ):

        chunk.rename(
            columns={chunk.columns[0]: "site_id"},
            inplace=True
        )

        total_input += len(chunk)

        keep_cols = ["site_id"] + [
            c for c in chunk.columns
            if c in valid_samples
        ]

        chunk = chunk[keep_cols]

        chunk = chunk[
            chunk["site_id"].isin(valid_sites)
        ].copy()

        total_after_manifest += len(chunk)

        if chunk.empty:
            print(
                f"Chunk {i + 1} | empty after filtering manifest"
            )
            continue

        before_missing_filter = len(chunk)

        chunk = filter_missing_values(chunk)

        after_missing_filter = len(chunk)
        removed_missing = before_missing_filter - after_missing_filter

        total_after_missing_filter += after_missing_filter

        if chunk.empty:
            print(
                f"Chunk {i + 1} | "
                f"0 CpGs after filtering missing values"
            )
            continue

        table = pa.Table.from_pandas(
            chunk,
            preserve_index=False
        )

        if writer is None:
            writer = pq.ParquetWriter(
                PARQUET_PATH,
                table.schema,
                compression="zstd"
            )

        writer.write_table(table)

        print(
            f"Chunk {i + 1} | "
            f"{after_missing_filter:,} CpGs saved | "
            f"{removed_missing:,} eliminated by >80% NA"
        )

    if writer:
        writer.close()

    print("\nFinal summary") 
    print("-------------") 
    print(f"CpGs read: {total_input:,}") 
    print(f"CpGs after manifest: {total_after_manifest:,}") 
    print(f"CpGs after NA filter: {total_after_missing_filter:,}") 
    print(f"File saved in: {PARQUET_PATH}")

# ============================================================
# Main
# ============================================================

if __name__ == "__main__":

    download_file()

    valid_sites, valid_samples = load_filters()

    process_methylation(
        valid_sites=valid_sites,
        valid_samples=valid_samples
    )
