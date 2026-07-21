# ============================================================
# Config
# ============================================================
import numpy as np
import pandas as pd
from pathlib import Path

# ============================================================
# Data Dir 
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "raw"
METHY_PATH = DATA_DIR / "methy_pb.parquet"
MANIFEST_PATH = DATA_DIR / "manifest_clean.parquet"
OUTPUT_PATH = DATA_DIR / "cpg_features.parquet"
MAX_NA_FRAC = 0.20
DTYPE = np.float32

# ============================================================
# Load Manifest
# ============================================================

print("Loading manifest...")
manifest = pd.read_parquet(
    MANIFEST_PATH,
    columns=["site_id"]
)

valid_sites = manifest["site_id"].values
print(
    f"Valid CpGs in manifest: "
    f"{len(valid_sites):,}"
)

# ============================================================
# Load Methylation
# ============================================================
print("Loading methylation matrix...")
df = pd.read_parquet(METHY_PATH)
print(df.shape)

# ============================================================
# Samples columns
# ============================================================
sample_cols = [
    c
    for c in df.columns
    if c != "site_id"
]

n_samples = len(sample_cols)

print(
    f"Samples: {n_samples}"
)

# ============================================================
# Filter Manifest
# ============================================================
print("Filtering manifest CpGs...")
df = df[
    df["site_id"].isin(valid_sites)
]
print(df.shape)

# ============================================================
# Matrix
# ============================================================
arr = df[sample_cols].to_numpy(
    dtype=DTYPE,
    copy=True
)

# ============================================================
# Filter >20% NA
# ============================================================

print("Filtering CpGs with >20% NA...")
na_frac = np.mean(
    np.isnan(arr),
    axis=1
)
valid_mask = (
    na_frac <= MAX_NA_FRAC
)
df = df.loc[
    valid_mask,
    ["site_id"]
].reset_index(drop=True)

arr = arr[valid_mask]
print(
    f"CpGs retained: "
    f"{arr.shape[0]:,}"
)

# ============================================================
# Imputation
# ============================================================

print("Imputing missing values...")
row_medians = np.nanmedian(
    arr,
    axis=1
)
nan_rows, nan_cols = np.where(
    np.isnan(arr)
)
arr[
    nan_rows,
    nan_cols
] = row_medians[
    nan_rows
]

# ============================================================
# Stats
# ============================================================

print("Computing statistics...")
median = np.median(
    arr,
    axis=1
).astype(DTYPE)

std = np.std(
    arr,
    axis=1,
    ddof=1
).astype(DTYPE)

# ============================================================
# Output
# ============================================================

cpg_feature = pd.DataFrame({
    "site_id": df["site_id"].values,
    "pb_median": median,
    "pb_std": std,
    "n_samples": n_samples
})

# ============================================================
# Save
# ============================================================
cpg_feature.to_parquet(
    OUTPUT_PATH,
    index=False,
    compression="zstd"
)

# ============================================================
# Summary
# ============================================================

print("SUMMARY:")
print(
    f"Final CpGs: "
    f"{len(cpg_feature):,}"
)
print(
    f"Samples: "
    f"{n_samples}"
)
print(
    f"Output: "
    f"{OUTPUT_PATH}"
)
print("\nDone.")
