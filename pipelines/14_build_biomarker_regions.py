"""
14_build_biomarker_regions_2.py
--------------------------------------------------

Build tumor-independent physical CpG regions from the
cleaned methylation manifest.

Purpose
-------
To generate a normalized region model for PostgreSQL and
Streamlit. Regions are defined from CpG genomic proximity and
manifest annotation, not from tumor-specific methylation,
expression, leukocyte methylation, or biological scores.

Inputs
------
- data/raw/manifest_clean.parquet

Processing
----------
1. Load and standardize manifest columns.
2. Normalize chromosomes and gene symbols.
3. Expand multi-gene CpG annotations when gene_mode='expand'.
4. Detect CpG clusters within each gene/chromosome group.
5. Keep clusters with at least --min-cpgs CpGs and adjacent
   CpGs separated by no more than --max-gap-bp.
6. Collapse duplicated physical intervals caused by multi-gene
   annotations.
7. Build a region table and a region-CpG bridge table.

Outputs
-------
- data/biomarkers/biomarker_region.parquet
- data/biomarkers/biomarker_region_cpg.parquet

Notes
-----
This step defines a reusable physical region universe. Tumor-
specific filters and biological prioritization are applied later.

--------------------------------------------------
"""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "raw"
OUT_DIR = ROOT / "data" / "biomarkers"


# ============================================================
# Helpers
# ============================================================


def first_existing_column(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    """Return the first existing column among candidate names."""
    for col in candidates:
        if col in df.columns:
            return col
    return None


def require_columns(df: pd.DataFrame, columns: Iterable[str], table_name: str) -> None:
    """Raise a helpful error if required columns are missing."""
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(
            f"Missing columns in {table_name}: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )


def normalize_chrom(chrom: Any) -> str:
    """Normalize chromosome labels by removing chr/CHR prefix."""
    if chrom is None:
        return ""
    try:
        if pd.isna(chrom):
            return ""
    except (TypeError, ValueError):
        pass
    value = str(chrom).strip()
    value = re.sub(r"^chr", "", value, flags=re.IGNORECASE)
    return value


def normalize_gene_symbol(gene: Any) -> str:
    """Normalize gene symbols for grouping."""
    if gene is None:
        return ""
    try:
        if pd.isna(gene):
            return ""
    except (TypeError, ValueError):
        pass
    value = str(gene).strip().upper()
    if value.lower() in {"", "nan", "none", "null", "na"}:
        return ""
    return value


def split_gene_symbols(value: Any) -> list[str]:
    """Split multi-gene annotations into unique normalized symbols."""
    if value is None:
        return []
    try:
        if pd.isna(value):
            return []
    except (TypeError, ValueError):
        pass

    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "null", "na"}:
        return []

    genes = [normalize_gene_symbol(part) for part in re.split(r"[;,]", text)]
    genes = [gene for gene in genes if gene]
    return list(dict.fromkeys(genes))


def is_low_priority_gene(gene_symbol: str) -> bool:
    """
    Heuristic used only to choose a display gene when several annotations map
    to the exact same physical region.

    The full list is still preserved in gene_symbols_all.
    """
    gene = normalize_gene_symbol(gene_symbol)
    low_priority_patterns = [
        r"^AC\d",       # e.g. AC009879.4
        r"^AL\d",
        r"^AP\d",
        r"^RP\d",
        r"^LOC\d",
        r"^LINC\d",
        r"^MIR\d",
        r"^SNOR",
        r"^RNU",
        r"^RNA5",
    ]
    return any(re.match(pattern, gene) for pattern in low_priority_patterns)


def sort_gene_symbols_for_display(genes: Iterable[Any]) -> list[str]:
    """Return unique genes ordered so protein-like symbols are shown first."""
    cleaned = [normalize_gene_symbol(g) for g in genes]
    cleaned = [g for g in cleaned if g]
    unique = list(dict.fromkeys(cleaned))
    return sorted(unique, key=lambda g: (is_low_priority_gene(g), g))


def join_gene_symbols(genes: Iterable[Any]) -> str:
    """Join unique gene symbols for compact annotation."""
    return ";".join(sort_gene_symbols_for_display(genes))


def make_region_id(chrom: str, start: int, end: int) -> str:
    """Create a stable short physical region ID from genomic coordinates only."""
    raw = f"{normalize_chrom(chrom)}|{int(start)}|{int(end)}"
    digest = hashlib.md5(raw.encode("utf-8")).hexdigest()[:10]
    return f"REG_{digest}"


def standardize_manifest_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize manifest columns to site_id, gene_symbol, chr, start_pos."""
    out = df.copy()
    rename: dict[str, str] = {}

    site_col = first_existing_column(out, ["site_id", "IlmnID", "illumina_id", "probe_id", "Name", "cpg_id"])
    gene_col = first_existing_column(out, ["gene_symbol", "gene", "Gene", "UCSC_RefGene_Name", "gene_name", "symbol"])
    chr_col = first_existing_column(out, ["chr", "chrom", "CHR", "chromosome", "Chromosome"])
    pos_col = first_existing_column(out, ["start_pos", "position", "pos", "MAPINFO", "start", "Start"])

    if site_col and site_col != "site_id":
        rename[site_col] = "site_id"
    if gene_col and gene_col != "gene_symbol":
        rename[gene_col] = "gene_symbol"
    if chr_col and chr_col != "chr":
        rename[chr_col] = "chr"
    if pos_col and pos_col != "start_pos":
        rename[pos_col] = "start_pos"

    if rename:
        out = out.rename(columns=rename)

    require_columns(out, ["site_id", "gene_symbol", "chr", "start_pos"], "manifest")
    return out


def prepare_manifest_rows(df: pd.DataFrame, gene_mode: str) -> pd.DataFrame:
    """Create one normalized row per CpG-gene-coordinate mapping."""
    out = standardize_manifest_columns(df)
    out = out[["site_id", "gene_symbol", "chr", "start_pos"]].copy()

    out["site_id"] = out["site_id"].astype(str)
    out["start_pos"] = pd.to_numeric(out["start_pos"], errors="coerce")
    out["chr"] = out["chr"].map(normalize_chrom)

    out = out.dropna(subset=["site_id", "gene_symbol", "chr", "start_pos"])
    out = out[(out["site_id"] != "") & (out["chr"] != "")].copy()
    out["start_pos"] = out["start_pos"].astype(int)

    input_rows = out.shape[0]

    if gene_mode == "first":
        out["gene_symbol"] = out["gene_symbol"].apply(lambda x: split_gene_symbols(x)[0] if split_gene_symbols(x) else "")
    else:
        out["gene_symbol"] = out["gene_symbol"].apply(split_gene_symbols)
        out = out.explode("gene_symbol").copy()

    out["gene_symbol"] = out["gene_symbol"].map(normalize_gene_symbol)
    out = out[out["gene_symbol"] != ""].copy()

    expanded_rows = out.shape[0]

    out = out.drop_duplicates(subset=["site_id", "gene_symbol", "chr", "start_pos"])

    print(f"Input manifest rows: {input_rows:,}")
    print(f"CpG-gene-coordinate rows after gene_mode='{gene_mode}': {expanded_rows:,}")
    print(f"Unique CpG-gene-coordinate rows used: {out.shape[0]:,}")
    print(f"Unique CpGs used: {out['site_id'].nunique():,}")
    print(f"Unique genes used: {out['gene_symbol'].nunique():,}")

    return out.reset_index(drop=True)


# ============================================================
# Region construction
# ============================================================


def finalize_gene_cluster(
    cluster: pd.DataFrame,
    flank_bp: int,
) -> tuple[dict[str, Any], pd.DataFrame]:
    """Convert one gene-specific CpG cluster into one candidate region row."""
    cluster = cluster.sort_values(["start_pos", "site_id"]).reset_index(drop=True)

    gene_symbol = str(cluster.loc[0, "gene_symbol"])
    chrom = str(cluster.loc[0, "chr"])
    core_start = int(cluster["start_pos"].min())
    core_end = int(cluster["start_pos"].max())
    core_length = int(core_end - core_start + 1)

    browser_start = max(1, int(core_start - flank_bp))
    browser_end = int(core_end + flank_bp)
    browser_length = int(browser_end - browser_start + 1)

    region_id = make_region_id(chrom, core_start, core_end)

    n_manifest_cpgs = int(cluster["site_id"].nunique())
    n_manifest_c_positions = int(cluster["start_pos"].nunique())
    cpg_density_per_100bp = float(n_manifest_cpgs / core_length * 100) if core_length > 0 else np.nan

    region_row = {
        "region_id": region_id,
        "gene_symbol": gene_symbol,
        "chr": chrom,
        "core_start": core_start,
        "core_end": core_end,
        "core_length": core_length,
        "browser_start": browser_start,
        "browser_end": browser_end,
        "browser_length": browser_length,
        "flank_bp": int(flank_bp),
        "n_manifest_cpgs": n_manifest_cpgs,
        "n_manifest_c_positions": n_manifest_c_positions,
        "cpg_density_per_100bp": cpg_density_per_100bp,
    }

    bridge = cluster[["site_id", "gene_symbol", "chr", "start_pos"]].copy()
    bridge.insert(0, "region_id", region_id)

    return region_row, bridge


def collapse_duplicate_physical_regions(
    candidate_regions: pd.DataFrame,
    candidate_region_cpg: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Collapse duplicated intervals caused by multi-gene annotations.

    The output region table has one row per physical interval. The primary
    display gene is stored in gene_symbol and all associated genes are kept in
    gene_symbols_all.
    """
    if candidate_regions.empty:
        empty_region = pd.DataFrame(columns=[
            "region_id", "gene_symbol", "gene_symbols_all", "n_associated_genes",
            "chr", "core_start", "core_end", "core_length",
            "browser_start", "browser_end", "browser_length", "flank_bp",
            "n_manifest_cpgs", "n_manifest_c_positions", "cpg_density_per_100bp",
            "region_rank_by_density",
        ])
        empty_bridge = pd.DataFrame(columns=[
            "region_id", "site_id", "gene_symbol", "site_gene_symbols_all", "chr", "start_pos", "cpg_order",
        ])
        return empty_region, empty_bridge

    key_cols = ["chr", "core_start", "core_end", "browser_start", "browser_end"]

    # Region-level physical table.
    region_rows: list[dict[str, Any]] = []
    for _, group in candidate_regions.groupby(key_cols, sort=False):
        first = group.iloc[0]
        genes = sort_gene_symbols_for_display(group["gene_symbol"].tolist())
        primary_gene = genes[0] if genes else ""
        gene_symbols_all = ";".join(genes)

        region_id = make_region_id(first["chr"], int(first["core_start"]), int(first["core_end"]))
        region_rows.append({
            "region_id": region_id,
            "gene_symbol": primary_gene,
            "gene_symbols_all": gene_symbols_all,
            "n_associated_genes": int(len(genes)),
            "chr": str(first["chr"]),
            "core_start": int(first["core_start"]),
            "core_end": int(first["core_end"]),
            "core_length": int(first["core_length"]),
            "browser_start": int(first["browser_start"]),
            "browser_end": int(first["browser_end"]),
            "browser_length": int(first["browser_length"]),
            "flank_bp": int(first["flank_bp"]),
        })

    regions = pd.DataFrame(region_rows)

    # Add region-level gene annotations to each candidate bridge row.
    region_gene_map = regions[[
        "region_id", "gene_symbol", "gene_symbols_all", "chr", "core_start", "core_end"
    ]].rename(columns={"gene_symbol": "region_gene_symbol", "gene_symbols_all": "region_gene_symbols_all"})

    bridge = candidate_region_cpg.merge(
        region_gene_map,
        on=["region_id", "chr"],
        how="left",
    )

    # Keep one row per physical region and CpG site. If the same site was linked
    # to the same physical interval through multiple gene annotations, preserve
    # those site-level annotations in site_gene_symbols_all.
    bridge_rows: list[dict[str, Any]] = []
    for (_, site_id), group in bridge.groupby(["region_id", "site_id"], sort=False):
        first = group.sort_values(["start_pos", "gene_symbol"]).iloc[0]
        site_genes = join_gene_symbols(group["gene_symbol"].tolist())
        bridge_rows.append({
            "region_id": first["region_id"],
            "site_id": str(site_id),
            "gene_symbol": first.get("region_gene_symbol", ""),
            "site_gene_symbols_all": site_genes,
            "chr": str(first["chr"]),
            "start_pos": int(first["start_pos"]),
        })

    region_cpg = pd.DataFrame(bridge_rows)
    region_cpg = region_cpg.sort_values(["region_id", "start_pos", "site_id"]).reset_index(drop=True)
    region_cpg["cpg_order"] = region_cpg.groupby("region_id").cumcount() + 1

    # Recalculate manifest counts after physical deduplication.
    counts = (
        region_cpg.groupby("region_id")
        .agg(
            n_manifest_cpgs=("site_id", "nunique"),
            n_manifest_c_positions=("start_pos", "nunique"),
        )
        .reset_index()
    )

    regions = regions.merge(counts, on="region_id", how="left")
    regions["n_manifest_cpgs"] = regions["n_manifest_cpgs"].fillna(0).astype(int)
    regions["n_manifest_c_positions"] = regions["n_manifest_c_positions"].fillna(0).astype(int)
    regions["cpg_density_per_100bp"] = np.where(
        regions["core_length"] > 0,
        regions["n_manifest_cpgs"] / regions["core_length"] * 100,
        np.nan,
    )

    regions = regions.sort_values(["chr", "core_start", "core_end", "gene_symbol"]).reset_index(drop=True)
    regions["region_rank_by_density"] = regions["cpg_density_per_100bp"].rank(
        method="dense", ascending=False
    ).astype(int)

    region_cpg = region_cpg.sort_values(["region_id", "cpg_order", "start_pos"]).reset_index(drop=True)

    return regions, region_cpg


def build_regions(
    manifest_rows: pd.DataFrame,
    max_gap_bp: int,
    min_cpgs: int,
    flank_bp: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build region and region-CpG tables from normalized manifest rows."""
    candidate_region_rows: list[dict[str, Any]] = []
    candidate_bridge_tables: list[pd.DataFrame] = []

    group_cols = ["gene_symbol", "chr"]

    for group_number, ((gene_symbol, chrom), group) in enumerate(manifest_rows.groupby(group_cols), start=1):
        group = group.sort_values(["start_pos", "site_id"]).reset_index(drop=True)
        if len(group) < min_cpgs:
            continue

        cluster_indices = [0]
        previous_pos = int(group.loc[0, "start_pos"])

        for idx in range(1, len(group)):
            current_pos = int(group.loc[idx, "start_pos"])
            gap = current_pos - previous_pos

            if gap <= max_gap_bp:
                cluster_indices.append(idx)
            else:
                if len(cluster_indices) >= min_cpgs:
                    region_row, bridge = finalize_gene_cluster(group.iloc[cluster_indices].copy(), flank_bp)
                    candidate_region_rows.append(region_row)
                    candidate_bridge_tables.append(bridge)
                cluster_indices = [idx]

            previous_pos = current_pos

        if len(cluster_indices) >= min_cpgs:
            region_row, bridge = finalize_gene_cluster(group.iloc[cluster_indices].copy(), flank_bp)
            candidate_region_rows.append(region_row)
            candidate_bridge_tables.append(bridge)

        if group_number % 5000 == 0:
            print(f"  Processed {group_number:,} gene/chromosome groups")

    if not candidate_region_rows:
        return collapse_duplicate_physical_regions(pd.DataFrame(), pd.DataFrame())

    candidate_regions = pd.DataFrame(candidate_region_rows)
    candidate_region_cpg = pd.concat(candidate_bridge_tables, ignore_index=True)

    n_candidate_regions = int(candidate_regions.shape[0])
    n_candidate_ids = int(candidate_regions["region_id"].nunique())

    regions, region_cpg = collapse_duplicate_physical_regions(candidate_regions, candidate_region_cpg)

    print(f"Candidate gene-specific regions before physical deduplication: {n_candidate_regions:,}")
    print(f"Unique physical region IDs before collapse: {n_candidate_ids:,}")
    print(f"Final physical regions after collapse: {regions.shape[0]:,}")
    print(f"Duplicated gene-annotated regions collapsed: {n_candidate_regions - regions.shape[0]:,}")

    return regions, region_cpg


# ============================================================
# CLI
# ============================================================


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build normalized physical CpG region tables from manifest coordinates."
    )

    parser.add_argument("--input-cpg", type=Path, default=RAW_DIR / "manifest_clean.parquet")
    parser.add_argument("--output-region", type=Path, default=OUT_DIR / "biomarker_region.parquet")
    parser.add_argument("--output-region-cpg", type=Path, default=OUT_DIR / "biomarker_region_cpg.parquet")

    parser.add_argument("--max-gap-bp", type=int, default=350)
    parser.add_argument("--min-cpgs", type=int, default=2)
    parser.add_argument("--flank-bp", type=int, default=100)
    parser.add_argument("--gene-mode", choices=["expand", "first"], default="expand")

    parser.add_argument("--write-csv", action="store_true", help="Also write CSV copies for inspection.")
    parser.add_argument("--limit", type=int, default=0, help="Optional input row limit for debugging.")

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.input_cpg.exists():
        raise FileNotFoundError(f"Missing input file: {args.input_cpg}")

    print(f"Loading manifest/CpG table: {args.input_cpg}")
    manifest = pd.read_parquet(args.input_cpg)

    if args.limit and args.limit > 0:
        manifest = manifest.head(args.limit).copy()
        print(f"Debug limit applied to input manifest: {manifest.shape[0]:,} rows")

    print(
        "Building tumor-independent physical regions with "
        f">={args.min_cpgs} CpGs, max adjacent gap {args.max_gap_bp} bp, "
        f"gene_mode='{args.gene_mode}'..."
    )

    manifest_rows = prepare_manifest_rows(manifest, gene_mode=args.gene_mode)

    regions, region_cpg = build_regions(
        manifest_rows=manifest_rows,
        max_gap_bp=args.max_gap_bp,
        min_cpgs=args.min_cpgs,
        flank_bp=args.flank_bp,
    )

    args.output_region.parent.mkdir(parents=True, exist_ok=True)
    args.output_region_cpg.parent.mkdir(parents=True, exist_ok=True)

    regions.to_parquet(args.output_region, index=False, compression="zstd")
    region_cpg.to_parquet(args.output_region_cpg, index=False, compression="zstd")

    print("\nSaved outputs:")
    print(f"  Region table:     {args.output_region}      shape={regions.shape}")
    print(f"  Region-CpG table: {args.output_region_cpg}  shape={region_cpg.shape}")

    if args.write_csv:
        region_csv = args.output_region.with_suffix(".csv")
        region_cpg_csv = args.output_region_cpg.with_suffix(".csv")
        regions.to_csv(region_csv, index=False)
        region_cpg.to_csv(region_cpg_csv, index=False)
        print(f"  Region CSV:       {region_csv}")
        print(f"  Region-CpG CSV:   {region_cpg_csv}")

    print("\nRegion columns:")
    for col in regions.columns:
        print(f"  - {col}")

    print("\nRegion-CpG columns:")
    for col in region_cpg.columns:
        print(f"  - {col}")


if __name__ == "__main__":
    main()
