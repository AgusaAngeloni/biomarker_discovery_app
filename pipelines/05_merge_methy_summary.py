"""
05_merge_methy_summary.py
--------------------------------------------------
Merges cohort-specific methylation summary tables
into a single pan-cancer dataset.

Processing includes:
    - Loading cohort summary files
    - Cohort annotation
    - Dataset concatenation
    - Variable standardization
    - Data type harmonization

Output:
    methylation_summary.parquet

This file represents the final analytical
dataset used for downstream analyses.
--------------------------------------------------
"""

# ============================================================
# Config
# ============================================================
from pathlib import Path
import pandas as pd

# ============================================================
# Data Dir and paths
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR  = ROOT / "data" 
INPUT_DIR = Path(DATA_DIR /"summary")

OUT_PATH =  DATA_DIR /"raw/methylation_summary.parquet"

# ============================================================
# Find Files
# ============================================================

files = list(INPUT_DIR.glob("*_summary.parquet"))
print(f"Found {len(files)} files")

# ============================================================
# Load + concat
# ============================================================

dfs = []
for file in files:
    print(f"Loading: {file.name}")
    # --------------------------------------------------------
    # Parse tumor type
    # Example:
    # COAD_summary.parquet
    # --------------------------------------------------------
    tumor_type = file.stem.split("_")[0]

    # --------------------------------------------------------
    # Load parquet
    # --------------------------------------------------------
    df = pd.read_parquet(file)

    # --------------------------------------------------------
    # Add tumor type
    # --------------------------------------------------------
    df["tumor_type"] = tumor_type
    dfs.append(df)

# ============================================================
# Concat
# ============================================================
final_df = pd.concat(
    dfs,
    ignore_index=True
)

# ============================================================
# Column Order
# ============================================================
final_df = final_df[
    [
        "site_id",     
        "tumor_type",     
        "tumor_median",   
        "tumor_std",
        "tumor_n", 
        "normal_median",
        "normal_std",   
        "normal_n",       
        "panTumor_median", 
        "panTumor_std",    
        "panTumor_n",      
        "panNormal_median",  
        "panNormal_std",    
        "panNormal_n",       
        "delta_median",     
        "HI_index",         
    ]
]

# ============================================================
# Types
# ============================================================
final_df = final_df.astype({
        "site_id": str ,
        "tumor_type": str,   
        "tumor_median":"float32",
        "tumor_std":"float32",
        "tumor_n":"int32",  
        "normal_median":"float32",
        "normal_std":"float32",
        "normal_n":"int32",  
        "panTumor_median":"float32",
        "panTumor_std":"float32",
        "panTumor_n":"int32",  
        "panNormal_median":"float32",
        "panNormal_std":"float32",
        "panNormal_n":"int32",  
        "delta_median":"float32",
        "HI_index": "float32"
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
