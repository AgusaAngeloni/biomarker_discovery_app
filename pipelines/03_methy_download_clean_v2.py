"""
03_methy_download_clean_v2.py
--------------------------------------------------
Downloads and filters TCGA HumanMethylation450K
beta-value data.

The methylation matrix is filtered using:
    - CpG probes present in the cleaned manifest
    - Samples present in the cleaned phenotype table

Processing is performed in chunks to minimize
memory consumption.

No probe-level quality filtering is applied at
this stage.

Output:
    methy.parquet
--------------------------------------------------
"""

# ============================================================
# Config
# ============================================================
from pathlib import Path
import requests
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# ============================================================
# Data Dir 
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "raw"
DATA_DIR.mkdir(parents=True, exist_ok=True)

METHYLATION_URL = (
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/"
    "GDC-PANCAN.methylation450.tsv"
)

TSV_PATH = DATA_DIR / "methylation_450k.tsv"
PARQUET_PATH = DATA_DIR / "methy.parquet"
MANIFEST_PATH = DATA_DIR / "manifest_clean.parquet"
PHENO_PATH = DATA_DIR / "pheno_clean.parquet"
CHUNK_SIZE = 5000

# ============================================================
# Functions 
# ============================================================
def download_file():
    if TSV_PATH.exists():
        print("✔ Archivo ya descargado")
        return

    print("Descargando archivo .tsv")
    r = requests.get(METHYLATION_URL, stream=True)
    r.raise_for_status()

    with open(TSV_PATH, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)

def load_filters():
    manifest = pd.read_parquet(MANIFEST_PATH)
    pheno = pd.read_parquet(PHENO_PATH)

    return (
        set(manifest["site_id"]),
        set(pheno["sample_id"])
    )

def process_methylation(valid_sites, valid_samples):
    writer = None
    for i, chunk in enumerate(
        pd.read_csv(TSV_PATH, sep="\t", chunksize=CHUNK_SIZE)
    ):
        chunk.rename(
            columns={chunk.columns[0]: "site_id"},
            inplace=True
        )
        keep_cols = ["site_id"] + [
            c for c in chunk.columns
            if c in valid_samples
        ]
        chunk = chunk[keep_cols]
        chunk = chunk[
            chunk["site_id"].isin(valid_sites)
        ]
        table = pa.Table.from_pandas(chunk)
        if writer is None:
            writer = pq.ParquetWriter(
                PARQUET_PATH,
                table.schema,
                compression="zstd"
            )

        writer.write_table(table)
        print(
            f"Chunk {i+1} | "
            f"{len(chunk):,} CpGs"
        )
    if writer:
        writer.close()

if __name__ == "__main__":
    download_file()
    valid_sites, valid_samples = load_filters()
    process_methylation(valid_sites, valid_samples)
