"""
06_build_gene_annotation.py
--------------------------------------------------
Downloads and preprocesses GENCODE human gene
annotation.

Purpose
-------
To generate a standardized gene-level annotation
table containing Ensembl Gene IDs, gene symbols,
genomic coordinates, strand information, biotype,
and transcription start site positions.

Inputs
------
- GENCODE v49 human annotation GTF

Processing
----------
1. Download the GENCODE GTF file if not available.
2. Load gene-level features only.
3. Select relevant genomic annotation fields.
4. Remove Ensembl version suffixes.
5. Calculate transcription start site coordinates
   according to gene strand.
6. Standardize column names.
7. Save the processed annotation as a Parquet file.

Outputs
-------
- gene_annotation.parquet

--------------------------------------------------
"""
#pip install pandas gtfparse pyarrow

# ============================================================
# Config
# ============================================================
from pathlib import Path
import requests
from gtfparse import read_gtf
import pandas as pd

# ============================================================
# Data Dir
# ============================================================
GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/"
    "Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
)
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR  = ROOT / "data" 

GTF_PATH = Path(DATA_DIR/"raw/gencode.v49.annotation.gtf.gz")
OUT_PATH = DATA_DIR/ "raw/gene_annotation.parquet"

# ============================================================
# Download
# ============================================================
GTF_PATH.parent.mkdir(parents=True, exist_ok=True)

if not GTF_PATH.exists():
    print("Downloading GTF...")
    response = requests.get(
        GTF_URL,
        stream=True
    )
    response.raise_for_status()

    with open(GTF_PATH, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    print("Download complete.")
else:
    print("GTF already exists.")

# ============================================================
# Load gene-level features only
# ============================================================
print("Loading GTF (genes only)...")
gtf = read_gtf(
    str(GTF_PATH),
    features={"gene"}
).to_pandas()

print("Loaded:", gtf.shape)

# ============================================================
# Select relevant genomic annotation fields
# ============================================================
genes = gtf[
    [
        "gene_id",
        "gene_name",
        "seqname",
        "start",
        "end",
        "strand",
        "gene_type",
    ]
]

# ============================================================
# Remove Ensembl version suffixes
# ============================================================
genes["gene_id"] = (
    genes["gene_id"]
    .str.replace(r"\..*", "")
)

# ============================================================
# Calculate transcription start site coordinates according to gene strand
# ============================================================
genes["tss"] = genes.apply(
    lambda x: x["start"]
    if x["strand"] == "+"
    else x["end"],
    axis=1
)

# ============================================================
# Standardize column names
# ============================================================
genes.columns = [
    "ensembl_id",
    "gene_symbol",
    "chr",
    "start_pos",
    "end_pos",
    "strand",
    "biotype",
    "tss",
]

# ============================================================
# Save
# ============================================================
genes.to_parquet(
    OUT_PATH,
    index=False
)

print("\nSaved:")
print(OUT_PATH)

print("\nPreview:")
print(genes.head())

print("\n✓ Processed successfully")
