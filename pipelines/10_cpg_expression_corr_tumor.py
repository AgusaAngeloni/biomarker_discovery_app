"""
10_cpg_expression_corr_tumor.py
--------------------------------------------------
Computes CpG-expression correlations in tumor
samples.

Purpose
-------
To evaluate the association between DNA
methylation and gene expression using
Spearman correlation.

Inputs
------
- methy.parquet
- expr.parquet
- pheno_clean.parquet
- cpg_gene_map.parquet

Processing
----------
1. Select tumor samples by cohort.
2. Match CpGs with associated genes.
3. Retrieve methylation beta-values.
4. Retrieve gene expression values.
5. Compute Spearman correlation.
6. Store correlation coefficients and p-values.

Outputs
-------
- correlations_tumor/{cohort}_Tumor_corr.parquet
--------------------------------------------------
"""
# ============================================================
# Config
# ============================================================
import gc
from itertools import product
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import stats

# ============================================================
# Data Dir
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
RAW_DIR         = ROOT / "data" / "raw"
METHY_PATH      = ROOT /RAW_DIR / "methy.parquet"
EXPR_PATH       = ROOT /RAW_DIR / "expr.parquet"
PHENO_PATH      = ROOT /RAW_DIR / "pheno_clean.parquet"
CPG_GENE_PATH   = ROOT /RAW_DIR / "cpg_gene_map.parquet"
OUTPUT_DIR      = Path(ROOT / "data" /"correlations_tumor")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

COHORTS         = ["COAD", "LUSC", "LIHC","LUAD"]
CLASSES         = ["Tumor"]
MIN_SAMPLES     = 10    # mínimo de samples para calcular correlación
DTYPE           = np.float32

# ============================================================
# Functions
# ============================================================
def strip_ensembl_version(eid: str) -> str:
    """ENSG00000123456.7 → ENSG00000123456"""
    return eid.split(".")[0]

def spearman_r(x: np.ndarray, y: np.ndarray):
    """Spearman r y p-value, ignorando NaN pareados."""
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() < MIN_SAMPLES:
        return np.nan, np.nan
    r, p = stats.spearmanr(x[mask], y[mask])
    return float(r), float(p)

# ============================================================
# Load
# ============================================================
print("Loading pheno...")
pheno = pd.read_parquet(PHENO_PATH).set_index("sample_id")
print(f"  {pheno.shape[0]:,} samples")

print("Loading cpg_gene_map...")
cpg_gene_map = pd.read_parquet(CPG_GENE_PATH)
print(f"  {len(cpg_gene_map):,} paires CpG↔Ensembl")

# Ensembl IDs we need (without version)
target_ensembl = set(cpg_gene_map["ensembl_id"].unique())

# ============================================================
# Methy parquet
# ============================================================
pf_methy      = pq.ParquetFile(METHY_PATH)
methy_cols    = set(pf_methy.schema_arrow.names)
n_rg          = pf_methy.metadata.num_row_groups
print(f"\nMethy: {pf_methy.metadata.num_rows:,} probes | {n_rg} row groups")

# ============================================================
# Expression parquet
# ============================================================
print("\nLoading expresión")
pf_expr    = pq.ParquetFile(EXPR_PATH)
expr_cols = pf_expr.schema_arrow.names  # parquet column names expr

# Detect gene column in expr
gene_col_expr = pf_expr.schema_arrow.names[0]  # assume first column = gene id

# Load the full expr — we filter genes and samples afterward.
# If it doesn't fit in RAM, see the note at the end.
expr = pf_expr.read().to_pandas().set_index(gene_col_expr)

# Strip version of the index: ENSG.X → ENSG
expr.index = expr.index.map(strip_ensembl_version)

# Filter only genes that are in cpg_gene_map
expr = expr.loc[expr.index.isin(target_ensembl)]
print(f"  Genes relevantes cargados: {len(expr):,}")
print(f"  Samples en expr: {expr.shape[1]:,}")

# Samples TCGA in expr (format TCGA-XX-XXXX-XXA)
tcga_mask = expr.columns.str.startswith("TCGA")
expr      = expr.loc[:, tcga_mask]
print(f"  Samples TCGA en expr: {expr.shape[1]:,}")

# Convert a float32
expr = expr.astype(DTYPE)
gc.collect()

# Main loop: per cohort × class
for cohort, sample_class in product(COHORTS, CLASSES):
    out_path = OUTPUT_DIR / f"{cohort}_{sample_class}_corr.parquet"
    if out_path.exists():
        print(f"\n⏭  {cohort} {sample_class} it already exists, jumping")
        continue

    print(f"\n{'='*55}")
    print(f"  {cohort} — {sample_class}")

    # Samples from the group in pheno
    mask = (
        (pheno["tissue_type"]  == cohort) &
        (pheno["sample_class"] == sample_class)
    )
    group_samples = pheno.index[mask].tolist()
    group_samples = [s for s in group_samples if s in methy_cols]
    print(f"  Samples in methy: {len(group_samples)}")

    if len(group_samples) < MIN_SAMPLES:
        print(f"  ⚠ Less than {MIN_SAMPLES} samples, skipping.")
        continue

    # Common samples between methy and expr for this group
    expr_sample_cols = set(expr.columns)
    common_samples   = [s for s in group_samples if s in expr_sample_cols]
    print(f"  Common samples methy∩expr: {len(common_samples)}")

    if len(common_samples) < MIN_SAMPLES:
        print(f"  ⚠ Less than {MIN_SAMPLES} common samples, skipping.")
        continue

    # Expression submatrix for these samples
    expr_sub = expr[common_samples]   # (n_genes, n_samples)

    # CpGs relevant to this run
    cpg_map_sub = cpg_gene_map.copy()

    # ── Read methylation by row groups and accumulate ─────────────────────────────
    # We build a dict: site_id → array of β-values ​​(only common_samples)

    print("  Reading methy by row groups")
    methy_dict = {}   # site_id → np.array float32 shape (n_common,)

    needed_sites = set(cpg_map_sub["site_id"].unique())
    needed_cols  = ["site_id"] + common_samples

    for rg_idx in range(n_rg):
        if rg_idx % 10 == 0:
            print(f"    Row group {rg_idx+1}/{n_rg}...", end="\r")

        table = pf_methy.read_row_group(rg_idx, columns=needed_cols)
        df_rg = table.to_pandas()
        del table

        # Filter only the CpGs we need
        df_rg = df_rg[df_rg["site_id"].isin(needed_sites)]

        for _, row in df_rg.iterrows():
            sid = row["site_id"]
            methy_dict[sid] = row[common_samples].values.astype(DTYPE)

        del df_rg
        gc.collect()

    print(f"\n  Loaded CpGs: {len(methy_dict):,}")

    # Calculate correlations 
    print("  Calculating Spearman correlations...")
    results = []

    for ensembl_id, group in cpg_map_sub.groupby("ensembl_id"):
        # Expression vector for this gene
        if ensembl_id not in expr_sub.index:
            continue
        expr_vec = expr_sub.loc[ensembl_id].values.astype(np.float64)

        for _, row in group.iterrows():
            site_id = row["site_id"]
            if site_id not in methy_dict:
                continue
            methy_vec = methy_dict[site_id].astype(np.float64)
            r, p      = spearman_r(methy_vec, expr_vec)
            results.append({
                "site_id":     site_id,
                "ensembl_id":  ensembl_id,
                "gene_symbol": row["gene_symbol"],
                "spearman_r":  r,
                "pvalue":      p,
                "n_samples":   int((~np.isnan(methy_vec) & ~np.isnan(expr_vec)).sum()),
            })
    df_result = pd.DataFrame(results)
    print(f"  Correlated pairs: {len(df_result):,}")

    if len(df_result) > 0:
        df_result.to_parquet(out_path, index=False, compression="zstd")
        print(f"  ✓ Saved: {out_path}")
        print(f"  Preview:\n{df_result.head(4).to_string()}\n")

    del methy_dict, expr_sub, df_result, results
    gc.collect()

print("\nDone. Archives in:", OUTPUT_DIR)

# ── NOTE on Expression RAM ───────────────────────────────────────────────
# If expr.parquet doesn't fit in RAM after filtering genes:
# Replace the entire load with:
#
# genes_needed = list(target_ensembl)
# expr = pf_expr.read(columns=["gene_id"] + samples_TCGA).to_pandas()...
#
# OR process genes in blocks with iter_batches().
