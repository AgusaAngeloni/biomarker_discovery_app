"""
04_methy_summary_v2.py
--------------------------------------------------
Computes cohort-level methylation summary statistics.

For each CpG site, the script calculates:
    - Median beta value
    - Standard deviation
    - Number of valid observations

Statistics are generated for:
    - Tumor samples
    - Normal samples
    - Pan-cancer tumor reference
    - Pan-cancer normal reference

A minimum completeness threshold is applied
before summary statistics are retained.

Output:
    {COHORT}_summary.parquet
--------------------------------------------------
"""
# ============================================================
# Config
# ============================================================
import gc
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

# ============================================================
# Data Dir
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "raw"

METHY_PATH = DATA_DIR / "methy.parquet"
PHENO_PATH = DATA_DIR / "pheno_clean.parquet"

OUTPUT_DIR = ROOT / "data" / "summary"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DTYPE = np.float32
MIN_VALID_FRACTION = 0.20

COHORTS = ["COAD", "LUSC", "LUAD", "LIHC"]
# -- Options for all tumor type
#COHORTS = sorted(
#    pheno["tissue_type"].unique()
#)

# ============================================================
# Functions 
# ============================================================
def load_pheno():
    return pd.read_parquet(PHENO_PATH).set_index("sample_id")

def classify_samples(pheno):
    tumor_mask = pheno["sample_class"] == "Tumor"
    normal_mask = pheno["sample_class"] == "Normal"
    groups = {}
    for cohort in COHORTS:
        mask = pheno["tissue_type"] == cohort
        groups[cohort] = {
            "tumor": pheno.index[mask & tumor_mask].tolist(),
            "normal": pheno.index[mask & normal_mask].tolist(),
        }
    groups["_all"] = {
        "tumor": pheno.index[tumor_mask].tolist(),
        "normal": pheno.index[normal_mask].tolist(),
    }
    return groups

def compute_stats(arr):
    n_samples = arr.shape[1]
    if n_samples == 0:
        nan = np.full(arr.shape[0], np.nan, dtype=np.float32)
        zero = np.zeros(arr.shape[0], dtype=np.int32)
        return nan, nan, zero
    n_valid = np.sum(~np.isnan(arr), axis=1)
    min_valid = max(
        1,
        int(n_samples * MIN_VALID_FRACTION)
    )
    median = np.nanmedian(arr, axis=1).astype(np.float32)
    std = np.nanstd(arr, axis=1, ddof=1).astype(np.float32)
    median[n_valid < min_valid] = np.nan
    std[n_valid < min_valid] = np.nan
    return median, std, n_valid.astype(np.int32)

def main():
    print(f"Loading phenotype")
    pheno = load_pheno()
    groups = classify_samples(pheno)
    print(f"Loading methylation")
    pf = pq.ParquetFile(METHY_PATH)
    methy_cols = set(pf.schema_arrow.names)

    for cohort in COHORTS:
        print(f"Procesing {cohort}:")
        cohort_tumor = set(groups[cohort]["tumor"])
        cohort_normal = set(groups[cohort]["normal"])
        group_cols = {
            "tumor": [s for s in groups[cohort]["tumor"] if s in methy_cols],
            "normal": [s for s in groups[cohort]["normal"] if s in methy_cols],
            "pan_tumor": [
                s for s in groups["_all"]["tumor"]
                if s in methy_cols and s not in cohort_tumor
            ],
            "pan_normal": [
                s for s in groups["_all"]["normal"]
                if s in methy_cols and s not in cohort_normal
            ]
        }
        all_needed = {"site_id"}
        for cols in group_cols.values():
            all_needed.update(cols)
        all_needed = [
            c for c in pf.schema_arrow.names
            if c in all_needed
        ]
        acc = {
            grp: {
                "median": [],
                "std": [],
                "_n": []
            }
            for grp in group_cols
        }
        probe_ids = []
        for rg in range(pf.metadata.num_row_groups):

            table = pf.read_row_group(
                rg,
                columns=all_needed
            )
            df = table.to_pandas()
            probe_ids.append(
                df["site_id"].values
            )
            for grp, cols in group_cols.items():
                arr = df[cols].values.astype(DTYPE)
                median, std, n_valid = compute_stats(arr)
                acc[grp]["median"].append(median)
                acc[grp]["std"].append(std)
                acc[grp]["_n"].append(n_valid)
            del df
            gc.collect()
        out = pd.DataFrame({
            "site_id": np.concatenate(probe_ids)
        })
        for grp in group_cols:
            out[f"{grp}_median"] = np.concatenate(
                acc[grp]["median"]
            )
            out[f"{grp}_std"] = np.concatenate(
                acc[grp]["std"]
            )
            out[f"{grp}_n"] = np.concatenate(
                acc[grp]["_n"]
            )
        out["delta_median"] = (
            out["tumor_median"] - out["normal_median"]
        ).astype(np.float32)
        denominator = np.sqrt(
            out["tumor_std"] ** 2 + out["normal_std"] ** 2
        )

        out["hi_index"] = np.where(
            denominator > 0,
            np.abs(out["delta_median"]) / denominator,
            np.nan
        ).astype(np.float32)
        out.to_parquet(
            OUTPUT_DIR / f"{cohort}_summary.parquet",
            index=False,
            compression="zstd"
        )
        print(
            f"✔ {cohort}: {out.shape}"
        )

if __name__ == "__main__":
    main()
