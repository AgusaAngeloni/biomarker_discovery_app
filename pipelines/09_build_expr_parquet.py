"""
09_build_expr_parquet.py
--------------------------------------------------
Downloads and preprocesses TCGA RNA-seq expression data.

Purpose
-------
To generate a filtered and compressed expression matrix containing 
only genes associated with CpGs retained in the methylation workflow.

Inputs
------
- TCGA Pan-Cancer HTSeq FPKM-UQ expression matrix
- pheno_clean.parquet
- cpg_gene_map.parquet

Processing
----------
1. Download TCGA expression data.
2. Retain valid TCGA samples.
3. Retain genes associated with CpGs.
4. Remove Ensembl version numbers.
5. Apply log2(FPKM-UQ + 1) transformation.
6. Store values as float32.

Outputs
-------
- expr.parquet
--------------------------------------------------
"""

# ============================================================
# Config
# ============================================================
import os
import gc
import requests
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path

# ============================================================
# Data Dir
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR  = ROOT / "data" / "raw"
EXPR_URL    = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/GDC-PANCAN.htseq_fpkm-uq.tsv"
TSV_PATH    = DATA_DIR / "expression_fpkm_uq.tsv"
OUTPUT_PATH = DATA_DIR / "expr.parquet"
CHUNK_SIZE  = 2000    # filas (genes) por chunk
DTYPE       = np.float32

# ============================================================
# Load
# ============================================================
print("Loading pheno...")
pheno = pd.read_parquet(DATA_DIR / "pheno_clean.parquet")
valid_samples = set(pheno["sample_id"])
print(f"  Valid samples: {len(valid_samples):,}")

print("Loading cpg_gene_map...")
cpg_gene_map = pd.read_parquet(DATA_DIR / "cpg_gene_map.parquet")

# The Ensembl files in cpg_gene_map are no longer versioned
valid_genes_no_version = set(cpg_gene_map["ensembl_id"].unique())
print(f"  Genes relevantes (sin versión): {len(valid_genes_no_version):,}")

if TSV_PATH.exists():
    print(f"\nTSV already exists in {TSV_PATH}, skipping download.")
else:
    print(f"\n⬇ Downloading expression (~4-8 GB, may take a while)...")
    response = requests.get(EXPR_URL, stream=True)
    response.raise_for_status()
    downloaded = 0
    with open(TSV_PATH, "wb") as f:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)
            downloaded += len(chunk)
            if downloaded % (100 * 1024 * 1024) == 0:
                print(f"  {downloaded // (1024*1024):,} MB downloaded", end="\r")
    print(f"\n✔ Saved in: {TSV_PATH}")

# ============================================================
# TSV Conversion → Parquet
# ============================================================
print("\n Converting TSV → Parquet with log2(FPKM-UQ + 1)\n")
writer      = None
total_rows  = 0
kept_rows   = 0

for i, chunk in enumerate(
    pd.read_csv(TSV_PATH, sep="\t", chunksize=CHUNK_SIZE)
):
    # Rename first column → gene_id
    chunk.rename(columns={chunk.columns[0]: "gene_id"}, inplace=True)

    # Strip version Ensembl: ENSG00000123.4 → ENSG00000123
    chunk["gene_id"] = chunk["gene_id"].str.split(".").str[0]

    # Filter genes: only those with associated CpGs
    chunk = chunk[chunk["gene_id"].isin(valid_genes_no_version)]

    if len(chunk) == 0:
        total_rows += CHUNK_SIZE
        continue

    # Filter columns: only valid samples of pheno + gene_id
    keep_cols = ["gene_id"] + [
        c for c in chunk.columns if c != "gene_id" and c in valid_samples
    ]
    chunk = chunk[keep_cols]

    # Apply log2(x + 1) to numeric columns
    sample_cols = [c for c in chunk.columns if c != "gene_id"]
    chunk[sample_cols] = np.log2(
        chunk[sample_cols].values.astype(np.float64) + 1
    ).astype(DTYPE)

    total_rows += CHUNK_SIZE
    kept_rows  += len(chunk)

    print(
        f"✅ Chunk {i+1} | "
        f"Saved genes: {len(chunk):,} | "
        f"Accumulated: {kept_rows:,}"
    )

    table = pa.Table.from_pandas(chunk, preserve_index=False)

    if writer is None:
        writer = pq.ParquetWriter(
            OUTPUT_PATH,
            table.schema,
            compression="zstd"
        )
        # Report columns on initialization
        n_samples_out = len(keep_cols) - 1
        print(f"\n📊 Samples in expr after filter: {n_samples_out:,}")
        print(f"   Writer parquet initialized\n")

    writer.write_table(table)
    del chunk, table
    gc.collect()

if writer:
    writer.close()

print(f"\n✓ Saved: {OUTPUT_PATH}")
print(f"  Total genes in processed TSVs: {total_rows:,}")
print(f"  Genes stored (with CpGs): {kept_rows:,}")

# ============================================================
# Quick verification
# ============================================================

print("\nChecking parquet...")
pf = pq.ParquetFile(OUTPUT_PATH)
print(f"  Rows (genes): {pf.metadata.num_rows:,}")
print(f"  Cols (gene_id + samples): {len(pf.schema_arrow.names):,}")
print(f"  Row groups: {pf.metadata.num_row_groups}")

sample_df = pf.read_row_group(0).to_pandas()
print(f"\n  Preview (first 3 rows, 5 columns):")
print(sample_df.iloc[:3, :5].to_string())
