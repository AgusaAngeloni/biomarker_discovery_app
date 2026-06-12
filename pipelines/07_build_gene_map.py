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
import pandas as pd
import requests
from pathlib import Path
import mygene
mg = mygene.MyGeneInfo()

# ============================================================
# Data Dir
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR  = ROOT / "data" / "raw"
MANIFEST_PATH = DATA_DIR/ "manifest_clean.parquet"
OUTPUT_PATH   = DATA_DIR/ "gene_map.parquet"
BATCH_SIZE    = 1000   # mygene.info acepta hasta 1000 por request
SPECIES       = "human"

# ============================================================
# Extract unique gene symbols from the manifest 
# ============================================================
print("Cargando manifest...")
manifest = pd.read_parquet(MANIFEST_PATH, columns=["gene"])

# ── Split por ';' ─────────────────────────────────────────────────────────────
symbols = (
    manifest["gene"]
    .dropna()
    .str.split(";")
    .explode()
    .str.strip()
    .loc[lambda s: s != ""]
    .unique()
    .tolist()
)
print(f"  Símbolos únicos en manifest: {len(symbols):,}")

# ── Query mygene.info in batches ──────────────────────────────────────────────

def query_mygene_batch(symbols_batch: list) -> list[dict]:
    """
    Consult mygene.info for a batch of gene symbols.
    """
    url = "https://mygene.info/v3/querymany"
    body = {
        "queries": symbols_batch,
        "scopes": "symbol",
        "fields": "symbol,ensembl.gene,entrezgene,type_of_gene",
        "species": "human",
        "size": 1,
        "dotfield": True,
    }
    resp = requests.post(url, json=body, timeout=60)
    resp.raise_for_status()
    return resp.json()

records = []
n_batches = (len(symbols) + BATCH_SIZE - 1) // BATCH_SIZE

print(f"\nConsulting mygene.info ({n_batches} batches of {BATCH_SIZE})...")

for i in range(n_batches):
    batch = symbols[i * BATCH_SIZE : (i + 1) * BATCH_SIZE]
    print(f"  Batch {i+1}/{n_batches}...", end="\r")

    try:
        results = mg.querymany(
            batch,
            scopes="symbol",
            fields="ensembl.gene,symbol",
            species="human",
            as_dataframe=False
        )
    except Exception as e:
        print(f"\n  ⚠ Error in batch {i+1}: {e}")
        time.sleep(5)

    results = mg.querymany(
        batch,
        scopes="symbol",
        fields="ensembl.gene,symbol",
        species="human",
        as_dataframe=False
    )

    for hit in results:
        # Ignore hits without match or Ensembl
        if hit.get("notfound") or "ensembl" not in hit:
            continue
        symbol = hit.get("symbol", hit.get("query", ""))
        ensembl_data = hit.get("ensembl")
        if not ensembl_data:
            continue

        # Like dict or list
        if isinstance(ensembl_data, dict):
            ensembl_data = [ensembl_data]
        ensembl_list = []
        for item in ensembl_data:
            gene_id = item.get("gene")
            if gene_id:
                ensembl_list.append(gene_id)
        for eid in ensembl_list:
            records.append({
                "gene_symbol": symbol,
                "ensembl_id":  eid,   # without version: ENSG00000XXXXXX
            })
    time.sleep(0.1)  # rate limit

print(f"\n  Hits with Ensembl: {len(records):,}")

# ============================================================
# Save
# ============================================================
gene_map = (
    pd.DataFrame(records)
    .drop_duplicates()
    .reset_index(drop=True)
)

print(f"  Unique pared symbol↔ensembl: {len(gene_map):,}")
print(f"  Preview:\n{gene_map.head(5).to_string()}\n")

gene_map.to_parquet(OUTPUT_PATH, index=False, compression="zstd")
print(f"✓ Saved: {OUTPUT_PATH}")
