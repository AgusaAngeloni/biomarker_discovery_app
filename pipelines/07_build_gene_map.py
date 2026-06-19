"""
07_build_gene_map.py
--------------------------------------------------
Builds a gene identifier conversion table linking
gene symbols to Ensembl Gene IDs.

Purpose
-------
To generate a standardized mapping between gene
symbols present in the methylation annotation
manifest and Ensembl identifiers used in TCGA
expression datasets.

Inputs
- manifest_clean.parquet

Processing
----------
1. Extract unique gene symbols from the manifest.
2. Query the MyGene.info API.
3. Retrieve corresponding Ensembl Gene IDs.
4. Remove duplicated mappings.
5. Generate a harmonized gene reference table.

Outputs
-------
- gene_map.parquet
--------------------------------------------------
"""
# ============================================================
# Config
# ============================================================
import time
from pathlib import Path
import pandas as pd
import mygene

# ============================================================
# Data Dir
# ============================================================

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "raw"

MANIFEST_PATH = DATA_DIR / "manifest_clean.parquet"
OUTPUT_PATH = DATA_DIR / "gene_map.parquet"

BATCH_SIZE = 1000
SPECIES = "human"

mg = mygene.MyGeneInfo()

# ============================================================
# Load manifest genes
# ============================================================

print("Loading manifest...")

manifest = pd.read_parquet(
    MANIFEST_PATH,
    columns=["gene"]
)

symbols = (
    manifest["gene"]
    .dropna()
    .astype(str)
    .str.split(";")
    .explode()
    .str.strip()
    .loc[lambda s: s != ""]
    .drop_duplicates()
    .tolist()
)

print(f"Unique gene symbols in manifest: {len(symbols):,}")

# ============================================================
# Query MyGene.info
# ============================================================

records = []

n_batches = (
    len(symbols) + BATCH_SIZE - 1
) // BATCH_SIZE

print(
    f"\nQuerying MyGene.info: "
    f"{n_batches} batches of {BATCH_SIZE}"
)

for i in range(n_batches):

    batch = symbols[
        i * BATCH_SIZE : (i + 1) * BATCH_SIZE
    ]

    print(
        f"Batch {i + 1}/{n_batches}",
        end="\r"
    )

    try:
        results = mg.querymany(
            batch,
            scopes="symbol,alias,prev_symbol",
            fields="symbol,ensembl.gene,entrezgene,type_of_gene",
            species=SPECIES,
            as_dataframe=False,
            returnall=False,
            verbose=False
        )

    except Exception as e:
        print(f"\nWARNING: error in batch {i + 1}: {e}")
        time.sleep(5)
        continue

    for hit in results:

        if hit.get("notfound"):
            continue

        # Original query symbol from manifest.
        # This is the symbol we must preserve for later merge.
        query_symbol = hit.get("query")

        # Current/official symbol returned by MyGene.info.
        official_symbol = hit.get(
            "symbol",
            query_symbol
        )

        ensembl_data = hit.get("ensembl")

        if not query_symbol or not ensembl_data:
            continue

        if isinstance(ensembl_data, dict):
            ensembl_data = [ensembl_data]

        for item in ensembl_data:

            ensembl_id = item.get("gene")

            if not ensembl_id:
                continue

            records.append({
                "gene_symbol": query_symbol,
                "official_symbol": official_symbol,
                "ensembl_id": ensembl_id
            })

    time.sleep(0.1)


# ============================================================
# Build output table
# ============================================================

gene_map = (
    pd.DataFrame(records)
    .drop_duplicates()
    .reset_index(drop=True)
)

print(f"\nHits with Ensembl: {len(gene_map):,}")

print(
    "Unique manifest symbols mapped: "
    f"{gene_map['gene_symbol'].nunique():,}"
)

print(
    "Unique Ensembl IDs: "
    f"{gene_map['ensembl_id'].nunique():,}"
)

print("\nPreview:")
print(gene_map.head(10).to_string(index=False))


# ============================================================
# Save
# ============================================================

gene_map.to_parquet(
    OUTPUT_PATH,
    index=False,
    compression="zstd"
)

print(f"\nSaved: {OUTPUT_PATH}")
