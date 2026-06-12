"""
13_generate_tumor_types.py
--------------------------------------------------
Generates a tumor-type reference table.

Purpose
-------
To create a lookup table containing TCGA
tumor abbreviations, full tumor names, and
tissue-of-origin annotations used by the
web application and database.

Inputs
------
None

Processing
----------
1. Define TCGA tumor codes.
2. Assign full cancer names.
3. Assign tissue categories.
4. Export reference table.

Outputs
-------
- tumor_types.parquet
--------------------------------------------------
"""
# ============================================================
# Config
# ============================================================
from pathlib import Path
import pandas as pd

# ============================================================
# Data
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR  = ROOT / "data" / "raw"
out_path = Path(DATA_DIR/"tumor_types.parquet")
out_path.parent.mkdir(parents=True, exist_ok=True)

TUMOR_TYPES = [
    # Liver disease — focus of the project
    ("LIHC",  "Liver Hepatocellular Carcinoma",               "liver"),
    ("CHOL",  "Cholangiocarcinoma",                           "liver"),
    ("LICA",  "Liver Intrahepatic Cholangiocarcinoma",        "liver"),

    # Gastrointestinal
    ("COAD",  "Colon Adenocarcinoma",                         "colon"),
    ("READ",  "Rectum Adenocarcinoma",                        "rectum"),
    ("STAD",  "Stomach Adenocarcinoma",                       "stomach"),
    ("ESCA",  "Esophageal Carcinoma",                         "esophagus"),
    ("PAAD",  "Pancreatic Adenocarcinoma",                    "pancreas"),

    # Urological
    ("KIRC",  "Kidney Renal Clear Cell Carcinoma",            "kidney"),
    ("BLCA",  "Bladder Urothelial Carcinoma",                 "bladder"),

    # Pulmonary
    ("LUAD",  "Lung Adenocarcinoma",                          "lung"),
    ("LUSC",  "Lung Squamous Cell Carcinoma",                 "lung"),

    # Other solids tumors
    ("BRCA",  "Breast Invasive Carcinoma",                    "breast"),
    ("PRAD",  "Prostate Adenocarcinoma",                      "prostate"),
    ("THCA",  "Thyroid Carcinoma",                            "thyroid"),
    ("UCEC",  "Uterine Corpus Endometrial Carcinoma",         "uterus"),
    ("SKCM",  "Skin Cutaneous Melanoma",                      "skin"),
    ("GBM",   "Glioblastoma Multiforme",                      "brain"),
]

# ============================================================
# Build Datafrane
# ============================================================
df = pd.DataFrame(TUMOR_TYPES, columns=["tumor_type", "full_name", "tissue"])
print(f"Rows       : {len(df)}")
print(f"Columns    : {list(df.columns)}")
print(f"\n{df.to_string(index=False)}\n")

# ============================================================
# Save
# ============================================================
df.to_parquet(out_path, index=False)

print(f"✅ Saved → {out_path}")
