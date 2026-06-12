"""
01_clean_manifiesto.py
--------------------------------------------------
Downloads and preprocesses the Illumina HumanMethylation450K
annotation manifest.

The script generates a standardized CpG reference table
used throughout the methylation analysis workflow.

Processing includes:
    - Download of hg38 manifest annotation
    - Selection of relevant genomic features
    - Retention of CpG probes only
    - Removal of sex chromosome and mitochondrial probes
    - Removal of probes lacking genomic annotation
    - Standardization of coordinate fields

Output:
    manifest_clean.parquet

--------------------------------------------------
"""

import os
import pandas as pd
from pathlib import Path

# ============================================================
# Config
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR  = ROOT / "data" / "raw"

os.makedirs(DATA_DIR, exist_ok=True)

# ============================================================
# Urls
# ============================================================

manifest_url = (
    "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.gencode.v36.tsv.gz"
)

# =========================
# Download of hg38 manifest annotation
# =========================

print("\n📖 Downloading Manifest HM4M50k...")

manifest = pd.read_csv(
    manifest_url,
    sep="\t",
    compression="gzip"
)
print(manifest)
print(manifest.head(1))

# =========================
# Selection of relevant genomic features
# =========================
print(manifest.columns.tolist())
manifest = manifest[
    [
        "genesUniq",
        "CpG_chrm",
        "probeID",
        "CpG_beg",
        "CpG_end",
        "distToTSS",
        "CGI",
        "CGIposition"
    ]
]

# =========================
# Retention of CpG probes only
# =========================

manifest = manifest[
    manifest["probeID"].str.startswith("cg")
]

# =========================
# Removal of sex chromosome and mitochondrial probes
# =========================


manifest = manifest[
    ~manifest["CpG_chrm"].isin(
        ["chrX", "chrY", "chrM"]
    )
]

manifest = manifest[
    manifest["CpG_chrm"].notna()
]

manifest = manifest[
    manifest["genesUniq"].notna()
]
manifest["CpG_beg"] = (
    manifest["CpG_beg"]
    .astype(str)
    .str.replace(",", "", regex=False)
)

manifest["CpG_beg"] = pd.to_numeric(
    manifest["CpG_beg"],
    errors="coerce"
)

manifest["CpG_beg"] = manifest["CpG_beg"].astype("Int64")

manifest["CpG_end"] = (
    manifest["CpG_end"]
    .astype(str)
    .str.replace(",", "", regex=False)
)

manifest["CpG_end"] = pd.to_numeric(
    manifest["CpG_end"],
    errors="coerce"
)

manifest["CpG_end"] = manifest["CpG_end"].astype("Int64")

# =========================
# RESET INDEX
# =========================

manifest = manifest.reset_index(drop=True)

# =========================
# Standardization of coordinate fields
# =========================
print("\n📖 Nulls data")
print(manifest.isna().sum())

print(f'New length: ',manifest.shape)

manifest.columns = [
    "gene",
    "chr",
    "site_id",
    "start",
    "end",
    "distance_tss",
    "cgi",
    "cg_island"
]

manifest.rename(columns={
    "gene":"gene",
    "chr":"chr",
    "site_id":"site_id",
    "start":"start_pos",
    "end":"end_pos",
    "distance_tss":"distance_tss",
    "cgi":"cgi",
    "cg_island":"cg_island"
}, inplace=True)

print(manifest.head())

# =========================
# Save Result
# =========================
manifest.to_parquet(
    os.path.join(
       DATA_DIR,
       "manifest_clean.parquet"),
    index=False
)

print("✔ Manifest clean save")
