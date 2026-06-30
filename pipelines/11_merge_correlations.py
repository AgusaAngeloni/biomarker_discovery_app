"""
11_merge_correlations.py
--------------------------------------------------
Merges cohort-specific methylation-expression
correlation results.

Purpose
-------
To generate a unified pan-cancer table
containing correlation metrics across all
tumor cohorts.

Inputs
------
- correlations_tumor/*.parquet

Processing
----------
1. Load cohort-specific correlation files.
2. Extract tumor type metadata.
3. Concatenate datasets.
4. Harmonize column types.

Outputs
-------
- expression_correlation.parquet
--------------------------------------------------
"""
# ============================================================
# Config
# ============================================================
from pathlib import Path
import pandas as pd

# ============================================================
# Data Dir
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
INPUT_DIR = Path(ROOT / "data" / "correlations_tumor")
OUT_PATH  = ROOT / "data" / "raw" / "expression_correlation.parquet"

# ============================================================
# Download
# ============================================================
files = list(INPUT_DIR.glob("*_corr.parquet"))
print(f"Found {len(files)} files")

# ============================================================
# Function
# ============================================================

dfs = []

for file in files:
    print(f"Loading: {file.name}")
    # --------------------------------------------------------
    # Parse filename
    # Example:
    # COAD_Tumor_corr.parquet
    # --------------------------------------------------------
    parts = file.stem.split("_")
    tumor_type   = parts[0]
   #sample_class = parts[1]
    # --------------------------------------------------------
    # Load parquet
    # --------------------------------------------------------
    df = pd.read_parquet(file)
    # --------------------------------------------------------
    # Add metadata columns
    # --------------------------------------------------------
    df["tumor_type"]   = tumor_type
    #  df["sample_class"] = sample_class
    dfs.append(df)

# ============================================================
# Concat
# ============================================================
final_df = pd.concat(
    dfs,
    ignore_index=True
)

# ============================================================
# Columns Order
# ============================================================
final_df = final_df[
    [
        "site_id",
        "tumor_type",
        #"sample_class",
        #"ensembl_id",
        #"gene_symbol",
        "spearman_r",
        "pvalue",
        "n_samples",
    ]
]

# ============================================================
# Types
# ============================================================
final_df = final_df.astype({
    "site_id": "string",
    "tumor_type": "string",
    #"sample_class": "string",
    #"ensembl_id": "string",
    #"gene_symbol": "string",
    "spearman_r": "float32",
    "pvalue": "float32",
    "n_samples": "int32",
})

# ============================================================
# Save
# ============================================================
final_df.to_parquet(
    OUT_PATH,
    index=False,
    compression="zstd"
)

print("\nSaved:")
print(OUT_PATH)

print("\nShape:")
print(final_df.shape)

print("\nPreview:")
print(final_df.head())
