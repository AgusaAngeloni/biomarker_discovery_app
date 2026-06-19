"""
08_build_cpg_gene_map.py
--------------------------------------------------
Builds a CpG-to-gene mapping table.

Purpose
-------
To establish relationships between CpG probes,
gene symbols, and Ensembl Gene IDs for
integrative methylation-expression analyses.

Inputs
------
- manifest_clean.parquet
- gene_map.parquet

Processing
----------
1. Expand CpGs annotated to multiple genes.
2. Generate unique CpG-gene pairs.
3. Merge with Ensembl identifiers.
4. Remove duplicated mappings.

Outputs
-------
- cpg_gene_map.parquet

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
DATA_DIR  = ROOT / "data" / "raw"
MANIFEST_PATH    = DATA_DIR/"manifest_clean.parquet"
GENE_MAP_PATH    = DATA_DIR/"gene_map.parquet"
OUTPUT_PATH      = DATA_DIR/"cpg_gene_map.parquet"

# ============================================================
# Download
# ============================================================
print("Loading manifest...")
manifest = pd.read_parquet(
    MANIFEST_PATH,
    columns=["site_id", "gene"]
)
print(f"  {len(manifest):,} CpGs")

print("Loadin gene_map...")
gene_map = pd.read_parquet(GENE_MAP_PATH)
print(f"  {len(gene_map):,} pairs symbol↔ensembl")

# ============================================================
# Expand CpGs annotated to multiple genes
# ============================================================
print("\nExpand the genes..")
cpg_genes = (
    manifest
    .dropna(subset=["gene"])
    .assign(
        gene_symbol=lambda df: df["gene"].str.split(";")
    )
    .explode("gene_symbol")
    .assign(gene_symbol=lambda df: df["gene_symbol"].str.strip())
    .loc[lambda df: df["gene_symbol"] != ""]
    .drop(columns=["gene"])
    .drop_duplicates(subset=["site_id", "gene_symbol"])   # duplicated isoforms
    .reset_index(drop=True)
)

print(f"  Unique Pairs CpG↔symbol: {len(cpg_genes):,}")

# ============================================================
# Join with gene_map for Ensembl
# ============================================================
print("\nJoin with Ensembl IDs...")
cpg_gene_map = (
    cpg_genes
    .merge(gene_map, on="gene_symbol", how="inner")  # inner: only those with Ensembl
    .drop_duplicates()
    .reset_index(drop=True)
)

print(f"  Unique pairs CpG↔Ensembl: {len(cpg_gene_map):,}")
print(f"  CpGs with at least one Ensembl: {cpg_gene_map['site_id'].nunique():,}")
print(f"  Unique genes mapped: {cpg_gene_map['ensembl_id'].nunique():,}")

# ── Estadísticas útiles ───────────────────────────────────────────────────────

cpgs_per_gene = cpg_gene_map.groupby("ensembl_id")["site_id"].count()
print(f"\n  CpGs per gene — median: {cpgs_per_gene.median():.0f} | "
      f"max: {cpgs_per_gene.max()} | min: {cpgs_per_gene.min()}")

print(f"\n  Preview:\n{cpg_gene_map.head(8).to_string()}\n")

# ── Guardar ───────────────────────────────────────────────────────────────────

cpg_gene_map.to_parquet(OUTPUT_PATH, index=False, compression="zstd")
print(f"✓ Saved: {OUTPUT_PATH}")
