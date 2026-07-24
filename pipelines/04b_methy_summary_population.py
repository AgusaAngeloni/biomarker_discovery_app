"""
04b_methy_summary_population.py
--------------------------------------------------
Computes CpG-level methylation summaries for population sensitivity analyses.

The default modes reproduce the preprint-oriented comparison:
    - full
    - asian_excluded

By default, population filtering is applied to the selected tumor cohort while
the pan-cancer reference remains unchanged. Use --pan-reference-mode matched to
apply the same population filter to the pan-cancer reference.

Inputs:
    data/raw/methy.parquet
    data/raw/pheno_population.parquet

Outputs:
    data/summary_population/{COHORT}__{MODE}__pan-{REFERENCE}_summary.parquet
--------------------------------------------------
"""

from __future__ import annotations

import argparse
import gc
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "raw"
METHY_PATH = DATA_DIR / "methy.parquet"
PHENO_PATH = DATA_DIR / "pheno_population.parquet"
OUTPUT_DIR = ROOT / "data" / "summary_population"

DEFAULT_COHORTS = ["COAD", "LUSC", "LUAD", "LIHC"]
DEFAULT_MODES = ["full", "asian_excluded"]
DTYPE = np.float32
MIN_VALID_FRACTION = 0.20


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute population-stratified CpG methylation summaries."
    )
    parser.add_argument("--methy", type=Path, default=METHY_PATH)
    parser.add_argument("--phenotype", type=Path, default=PHENO_PATH)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument(
        "--cohorts",
        nargs="+",
        default=DEFAULT_COHORTS,
        help="TCGA cohort abbreviations, for example LUSC LUAD.",
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        default=DEFAULT_MODES,
        choices=["full", "asian_excluded", "reported_non_asian", "asian_only"],
    )
    parser.add_argument(
        "--pan-reference-mode",
        choices=["full", "matched"],
        default="full",
        help=(
            "full keeps the original pan-cancer reference for every comparison; "
            "matched applies the selected population filter to the pan-cancer reference."
        ),
    )
    parser.add_argument(
        "--min-valid-fraction",
        type=float,
        default=MIN_VALID_FRACTION,
    )
    return parser.parse_args()


def validate_inputs(args: argparse.Namespace) -> None:
    for path in [args.methy, args.phenotype]:
        if not path.exists():
            raise FileNotFoundError(path)
    if not 0 < args.min_valid_fraction <= 1:
        raise ValueError("--min-valid-fraction must be in (0, 1].")


def load_phenotype(path: Path) -> pd.DataFrame:
    pheno = pd.read_parquet(path)
    required = {
        "sample_id",
        "tissue_type",
        "sample_class",
        "race_reported",
        "is_asian",
    }
    missing = sorted(required.difference(pheno.columns))
    if missing:
        raise KeyError(
            "Population phenotype is missing columns: " + ", ".join(missing)
        )

    pheno = pheno.copy()
    pheno["sample_id"] = pheno["sample_id"].astype("string").str.strip()
    pheno["race_reported"] = pheno["race_reported"].fillna(False).astype(bool)
    pheno["is_asian"] = pheno["is_asian"].fillna(False).astype(bool)

    if pheno["sample_id"].duplicated().any():
        duplicates = pheno.loc[
            pheno["sample_id"].duplicated(keep=False), "sample_id"
        ].head(10)
        raise ValueError(f"Duplicate sample_id values in phenotype: {duplicates.tolist()}")

    return pheno.set_index("sample_id", drop=False)


def population_mask(pheno: pd.DataFrame, mode: str) -> pd.Series:
    if mode == "full":
        return pd.Series(True, index=pheno.index)
    if mode == "asian_excluded":
        # Unknown race remains included, matching a literal exclusion analysis.
        return ~pheno["is_asian"]
    if mode == "reported_non_asian":
        return pheno["race_reported"] & ~pheno["is_asian"]
    if mode == "asian_only":
        return pheno["is_asian"]
    raise ValueError(f"Unsupported population mode: {mode}")


def sample_ids(
    pheno: pd.DataFrame,
    mask: pd.Series,
    sample_class: str,
    cohort: str | None = None,
    exclude_cohort: str | None = None,
) -> list[str]:
    selected = mask & pheno["sample_class"].eq(sample_class)
    if cohort is not None:
        selected &= pheno["tissue_type"].eq(cohort)
    if exclude_cohort is not None:
        selected &= ~pheno["tissue_type"].eq(exclude_cohort)
    return pheno.index[selected].tolist()


def compute_stats(
    arr: np.ndarray,
    min_valid_fraction: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_sites = arr.shape[0]
    n_samples = arr.shape[1]
    if n_samples == 0:
        nan = np.full(n_sites, np.nan, dtype=np.float32)
        zero = np.zeros(n_sites, dtype=np.int32)
        return nan, nan, zero

    n_valid = np.sum(~np.isnan(arr), axis=1)
    min_valid = max(1, int(n_samples * min_valid_fraction))

    with np.errstate(all="ignore"):
        median = np.nanmedian(arr, axis=1).astype(np.float32)
        std = np.nanstd(arr, axis=1, ddof=1).astype(np.float32)

    median[n_valid < min_valid] = np.nan
    std[n_valid < min_valid] = np.nan
    return median, std, n_valid.astype(np.int32)


def build_groups(
    pheno: pd.DataFrame,
    cohort: str,
    mode: str,
    pan_reference_mode: str,
    methy_columns: set[str],
) -> dict[str, list[str]]:
    target_mask = population_mask(pheno, mode)
    pan_mask = (
        population_mask(pheno, mode)
        if pan_reference_mode == "matched"
        else population_mask(pheno, "full")
    )

    groups = {
        "tumor": sample_ids(pheno, target_mask, "Tumor", cohort=cohort),
        "normal": sample_ids(pheno, target_mask, "Normal", cohort=cohort),
        "pan_tumor": sample_ids(
            pheno,
            pan_mask,
            "Tumor",
            exclude_cohort=cohort,
        ),
        "pan_normal": sample_ids(
            pheno,
            pan_mask,
            "Normal",
            exclude_cohort=cohort,
        ),
    }

    return {
        group: [sample for sample in samples if sample in methy_columns]
        for group, samples in groups.items()
    }


def summarize_one(
    pf: pq.ParquetFile,
    group_columns: dict[str, list[str]],
    cohort: str,
    mode: str,
    pan_reference_mode: str,
    min_valid_fraction: float,
) -> pd.DataFrame:
    all_needed = {"site_id"}
    for columns in group_columns.values():
        all_needed.update(columns)
    ordered_columns = [
        column for column in pf.schema_arrow.names if column in all_needed
    ]

    accumulators = {
        group: {"median": [], "std": [], "n": []}
        for group in group_columns
    }
    site_ids: list[np.ndarray] = []

    for row_group in range(pf.metadata.num_row_groups):
        frame = pf.read_row_group(row_group, columns=ordered_columns).to_pandas()
        site_ids.append(frame["site_id"].astype(str).to_numpy())

        for group, columns in group_columns.items():
            if columns:
                array = frame[columns].to_numpy(dtype=DTYPE, copy=False)
            else:
                array = np.empty((len(frame), 0), dtype=DTYPE)

            median, std, n_valid = compute_stats(array, min_valid_fraction)
            accumulators[group]["median"].append(median)
            accumulators[group]["std"].append(std)
            accumulators[group]["n"].append(n_valid)

        del frame
        gc.collect()

    output = pd.DataFrame({"site_id": np.concatenate(site_ids)})
    output.insert(1, "tumor_type", cohort)
    output.insert(2, "population_mode", mode)
    output.insert(3, "pan_reference_mode", pan_reference_mode)

    for group in group_columns:
        output[f"{group}_median"] = np.concatenate(accumulators[group]["median"])
        output[f"{group}_std"] = np.concatenate(accumulators[group]["std"])
        output[f"{group}_n"] = np.concatenate(accumulators[group]["n"])
        output[f"{group}_total_samples"] = len(group_columns[group])

    output["delta_median"] = (
        output["tumor_median"] - output["normal_median"]
    ).astype(np.float32)

    denominator = np.sqrt(
        output["tumor_std"] ** 2 + output["normal_std"] ** 2
    )
    output["hi_index"] = np.where(
        denominator > 0,
        np.abs(output["delta_median"]) / denominator,
        np.nan,
    ).astype(np.float32)
    output["dispersion_index"] = output["hi_index"]

    return output


def main() -> None:
    args = parse_args()
    validate_inputs(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading phenotype: {args.phenotype}")
    pheno = load_phenotype(args.phenotype)
    available_cohorts = set(pheno["tissue_type"].dropna().astype(str))
    unknown_cohorts = sorted(set(args.cohorts).difference(available_cohorts))
    if unknown_cohorts:
        raise ValueError(f"Cohorts absent from phenotype: {unknown_cohorts}")

    print(f"Opening methylation matrix: {args.methy}")
    pf = pq.ParquetFile(args.methy)
    if "site_id" not in pf.schema_arrow.names:
        raise KeyError("methy.parquet must contain a site_id column.")
    methy_columns = set(pf.schema_arrow.names)

    for cohort in args.cohorts:
        for mode in args.modes:
            print(
                f"\nProcessing cohort={cohort}, mode={mode}, "
                f"pan_reference={args.pan_reference_mode}"
            )
            groups = build_groups(
                pheno=pheno,
                cohort=cohort,
                mode=mode,
                pan_reference_mode=args.pan_reference_mode,
                methy_columns=methy_columns,
            )

            for group, samples in groups.items():
                print(f"  {group:10s}: {len(samples):4d} methylation samples")

            if not groups["tumor"]:
                print("  WARNING: no tumor samples; statistics will be NA.")
            if not groups["normal"]:
                print("  WARNING: no normal samples; delta and HI will be NA.")

            output = summarize_one(
                pf=pf,
                group_columns=groups,
                cohort=cohort,
                mode=mode,
                pan_reference_mode=args.pan_reference_mode,
                min_valid_fraction=args.min_valid_fraction,
            )

            output_path = args.output_dir / (
                f"{cohort}__{mode}__pan-{args.pan_reference_mode}_summary.parquet"
            )
            output.to_parquet(output_path, index=False, compression="zstd")
            print(f"Saved: {output_path}  shape={output.shape}")


if __name__ == "__main__":
    main()
