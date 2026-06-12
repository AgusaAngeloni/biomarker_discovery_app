"""
02_clean_phenotype.py
--------------------------------------------------
Downloads and preprocesses TCGA phenotype metadata.

The resulting table provides harmonized sample
annotations used to classify methylation profiles
according to tissue type and disease status.

Processing includes:
    - Download of TCGA phenotype metadata
    - Selection of TCGA samples
    - Retention of primary tumor and normal tissues
    - Extraction of demographic variables
    - Generation of patient identifiers
    - Cohort and sample class annotation

Output:
    pheno_clean.parquet
--------------------------------------------------
"""

# ============================================================
# Config
# ============================================================
import os
import pandas as pd
import requests
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path

# ============================================================
# Data Dir
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR  = ROOT / "data" / "raw"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Urls y Paths
# ============================================================
PHENOTYPE_URL = (
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/GDC-PANCAN.basic_phenotype.tsv"
)
pheno_gz = DATA_DIR / "phenotype.tsv"

# ============================================================
# Functions
# ============================================================
def download_file(url, output_path):

    print(f"\n⬇ Descargando:\n{url}")

    response = requests.get(url, stream=True)
    response.raise_for_status()

    with open(output_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)

    print(f"✔ Guardado en: {output_path}")

# ============================================================
# Download of TCGA phenotype metadata
# ============================================================
print("\n📖 Descargando phenotype...")

if not pheno_gz.exists():

    print("Downloading phenotype...")
    download_file(PHENOTYPE_URL, pheno_gz)
    print("Download complete.")

else:
    print("Phenotype already exists.")

print("\n📖 Leyendo phenotype...")

pheno_df = pd.read_csv(
    pheno_gz,
    sep="\t"
)

print(pheno_df.shape)

# ============================================================
# Selection of TCGA samples
# ============================================================

pheno_df = pheno_df[
    pheno_df["program"] == "TCGA"
].copy()

print("\n✔ Solo TCGA:")
print(pheno_df.shape)

# ============================================================
# Retention of primary tumor and normal tissues:
#   -01A = Primary Tumor
#   -11A = Solid Tissue Normal
# ============================================================

pheno_df = pheno_df[
    pheno_df["sample"].str.endswith(("-01A", "-11A"))
].copy()

print("\n✔ Solo 01A y 11A:")
print(pheno_df.shape)

# ============================================================
# Extraction of demographic variables
# ============================================================

keep_cols = [
    "sample",
    "project_id",
    "sample_type",
    "sample_type_id",
    "Gender",
    "Age at Diagnosis in Years"
]

pheno_df = pheno_df[keep_cols].copy()

# ----------- Rename
pheno_df.rename(columns={
    "sample": "sample_id",
    "project_id": "tumor_type",
    "sample_type": "sample_type",
    "sample_type_id": "sample_type_id",
    "Gender": "sex",
    "Age at Diagnosis in Years": "age"
}, inplace=True)

# ============================================================
# Generation of patient identifiers
# ============================================================
pheno_df["patient_id"] = (
    pheno_df["sample_id"]
    .str.split("-")
    .str[:3]
    .str.join("-")
)

# ---------- Conversion age → integer
pheno_df["age"] = pd.to_numeric(
    pheno_df["age"],
    errors="coerce"
).astype("Int32")

# ============================================================
# Cohort and sample class annotation
# ============================================================
pheno_df["tissue_type"] = pheno_df[
    "tumor_type"
].str.replace("TCGA-", "", regex=False)

sample_class_map = {
    1: "Tumor",
    11: "Normal"
}

pheno_df["sample_class"] = pheno_df[
    "sample_type_id"
].map(sample_class_map)

# =========================
# Save Result
# =========================
pheno_df.to_parquet(
    os.path.join(
        DATA_DIR,
       "pheno_clean.parquet"),
    index=False
)

print("✔ Pheno limpio guardado")
