"""
build_postgres_region_model_v2.py
───────────────────────────────────────────────────────────────
Create the PostgreSQL schema and import derived .parquet tables.

This version keeps the same database connection strategy used in the
project's build_postgres.py:

1) Local PostgreSQL using separate environment variables:

    export POSTGRES_USER=postgres
    export POSTGRES_PASSWORD=your_password
    export POSTGRES_HOST=localhost
    export POSTGRES_PORT=5432
    export POSTGRES_DB=db_methylation

    python pipelines/build_postgres_region_model.py

2) Direct database URL, useful for remote PostgreSQL providers:

    export DATABASE_URL="postgresql+psycopg2://user:password@host:5432/dbname"
    python pipelines/build_postgres_region_model.py

Notes:
    - If DATABASE_URL is provided, the script connects directly to that DB.
      It does not try to create the database itself.
    - If DATABASE_URL is not provided, the script builds the URL from POSTGRES_* variables
      and attempts to create POSTGRES_DB if it does not exist.
    - Existing project tables are dropped and recreated by default.
    - Large input files should exist under data/raw/ and data/biomarkers/.

Region model tables added/updated for the current GitHub app:
    - biomarker_cpg_score
    - biomarker_region
    - biomarker_region_cpg
    - biomarker_region_sequence_score

Compatibility notes:
    - Preserves gene_symbols_all and site_gene_symbols_all from pipeline 14.
    - Keeps both n_manifest_c and n_manifest_c_positions so existing pages work.
    - Preserves sequence_start/sequence_end/sequence_error and sequence density fields from pipeline 15.
    - Preserves gene-specific expression-correlation columns when available.

Requirements:
    pip install pandas pyarrow sqlalchemy psycopg2-binary
"""

from __future__ import annotations

from io import StringIO
from pathlib import Path
import argparse
import os
import re
import traceback
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd
from sqlalchemy import create_engine, text
from sqlalchemy.engine import URL


# ============================================================
# PATHS
# ============================================================

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "raw"
BIOMARKER_DIR = ROOT / "data" / "biomarkers"


# ============================================================
# DATABASE CONFIG — SAME LOGIC AS build_postgres.py
# ============================================================

DATABASE_URL = os.getenv("DATABASE_URL")

POSTGRES_USER = os.getenv("POSTGRES_USER", "postgres")
POSTGRES_PASSWORD = os.getenv("POSTGRES_PASSWORD", "")
POSTGRES_HOST = os.getenv("POSTGRES_HOST", "localhost")
POSTGRES_PORT = int(os.getenv("POSTGRES_PORT", "5432"))
POSTGRES_DB = os.getenv("POSTGRES_DB", "db_methylation")


def build_url(database: str) -> URL:
    """Build a SQLAlchemy PostgreSQL URL safely, including special characters in passwords."""
    return URL.create(
        drivername="postgresql+psycopg2",
        username=POSTGRES_USER,
        password=POSTGRES_PASSWORD,
        host=POSTGRES_HOST,
        port=POSTGRES_PORT,
        database=database,
    )


def get_engine():
    """
    Return a SQLAlchemy engine.

    This intentionally keeps the same connection strategy as build_postgres.py:
    1) If DATABASE_URL is present, use it directly.
    2) Otherwise, use POSTGRES_* variables, create POSTGRES_DB if needed,
       and connect to that database.
    """
    if DATABASE_URL:
        print("\n🔗 Using DATABASE_URL")
        print("   If you wanted to use POSTGRES_USER/POSTGRES_PASSWORD instead, run: unset DATABASE_URL")
        return create_engine(DATABASE_URL)

    print("\n🔍 Checking local PostgreSQL database using POSTGRES_* variables...\n")
    print(f"   user={POSTGRES_USER}")
    print(f"   host={POSTGRES_HOST}")
    print(f"   port={POSTGRES_PORT}")
    print(f"   db={POSTGRES_DB}")
    print(f"   password_set={bool(POSTGRES_PASSWORD)}")

    admin_engine = create_engine(build_url("postgres"))

    with admin_engine.connect() as conn:
        conn = conn.execution_options(isolation_level="AUTOCOMMIT")

        exists = conn.execute(
            text("SELECT 1 FROM pg_database WHERE datname = :dbname"),
            {"dbname": POSTGRES_DB},
        ).scalar()

        if not exists:
            print(f"   Creating database: {POSTGRES_DB}")
            conn.execute(text(f'CREATE DATABASE "{POSTGRES_DB}"'))
            print("   ✓ Database created")
        else:
            print(f"   ✓ Database already exists: {POSTGRES_DB}")

    admin_engine.dispose()
    return create_engine(build_url(POSTGRES_DB))


# ============================================================
# PARQUET INPUTS
# ============================================================

# Values may be a single Path or a list of fallback paths. The first existing
# file is loaded. This keeps compatibility with the previous and current names.
FILES: Dict[str, Path | List[Path]] = {
    # Base project tables
    "cpg_annotation": DATA_DIR / "manifest_clean.parquet",
    "cpg_gene_map": DATA_DIR / "cpg_gene_map.parquet",
    "expression_correlation": DATA_DIR / "expression_correlation.parquet",
    "cpg_features": DATA_DIR / "cpg_features.parquet",
    "sample_metadata": DATA_DIR / "pheno_clean.parquet",
    "gene_annotation": DATA_DIR / "gene_annotation.parquet",
    "tumor_types": DATA_DIR / "tumor_types.parquet",
    "tumor_summary": DATA_DIR / "methylation_summary.parquet",

    # Region biomarker model tables
    "biomarker_cpg_score": [
        BIOMARKER_DIR / "biomarker_cpg_score.parquet",
        BIOMARKER_DIR / "biomarker_cpg_biological_evidence.parquet",
    ],
    "biomarker_region": [
        BIOMARKER_DIR / "biomarker_region.parquet",
        BIOMARKER_DIR / "biomarker_region_universe.parquet",
    ],
    "biomarker_region_cpg": [
        BIOMARKER_DIR / "biomarker_region_cpg.parquet",
        BIOMARKER_DIR / "biomarker_region_universe.parquet",
    ],
    "biomarker_region_sequence_score": [
        BIOMARKER_DIR / "biomarker_region_sequence_score.parquet",
        BIOMARKER_DIR / "biomarker_region_sequence_features.parquet",
    ],
}

BASE_TABLES = [
    "cpg_annotation",
    "cpg_gene_map",
    "expression_correlation",
    "cpg_features",
    "sample_metadata",
    "gene_annotation",
    "tumor_types",
    "tumor_summary",
]

BIOMARKER_TABLES = [
    "biomarker_cpg_score",
    "biomarker_region",
    "biomarker_region_cpg",
    "biomarker_region_sequence_score",
]


# ============================================================
# TABLE COLUMN ORDER
# ============================================================

TABLE_COLUMNS: Dict[str, List[str]] = {
    "cpg_annotation": [
        "site_id",
        "gene",
        "chr",
        "start_pos",
        "end_pos",
        "distance_tss",
        "cgi",
        "cg_island",
        "relation_to_island",
        "ucsc_refgene_group",
    ],
    "cpg_gene_map": [
        "site_id",
        "ensembl_id",
        "gene_symbol",
    ],
    "expression_correlation": [
        "site_id",
        "tumor_type",
        "sample_class",
        "ensembl_id",
        "gene_symbol",
        "spearman_r",
        "pvalue",
        "n_samples",
    ],
    "cpg_features": [
        "site_id",
        "leukocyte_median",
        "leukocyte_std",
        "n_samples",
    ],
    "sample_metadata": [
        "sample_id",
        "patient_id",
        "tumor_type",
        "sample_type",
        "sample_type_id",
        "sex",
        "age",
        "tissue_type",
        "sample_class",
        "platform",
        "batch",
    ],
    "gene_annotation": [
        "ensembl_id",
        "gene_symbol",
        "chr",
        "start_pos",
        "end_pos",
        "strand",
        "biotype",
        "tss",
    ],
    "tumor_types": [
        "tumor_type",
        "full_name",
        "tissue",
    ],
    # This table keeps original column names used by older app pages and also
    # adds aliases used by the dynamic region scoring page.
    "tumor_summary": [
        "site_id",
        "tumor_type",
        "tumor_median",
        "tumor_std",
        "tumor_n",
        "n_tumor",
        "normal_median",
        "normal_std",
        "normal_n",
        "n_normal",
        "pan_tumor_median",
        "pan_tumor_std",
        "pan_tumor_n",
        "pan_normal_median",
        "pan_normal_std",
        "pan_normal_n",
        "delta_median",
        "hi_index",
        "dispersion_index",
    ],
    "biomarker_cpg_score": [
        "site_id",
        "tumor_type",
        "gene_symbol",
        "biological_score",
        "delta_score",
        "normal_low_score",
        "leukocyte_low_score",
        "pancancer_specificity_score",
        "hi_score",
        "expression_score",
        "sample_support_score",
        "annotation_score",
        "passes_loose_seed",
        "passes_default_filter",
        "passes_strict_filter",
    ],
    "biomarker_region": [
        "region_id",
        "gene_symbol",
        "gene_symbols_all",
        "n_associated_genes",
        "chr",
        "core_start",
        "core_end",
        "core_length",
        "browser_start",
        "browser_end",
        "browser_length",
        "flank_bp",
        "n_manifest_cpgs",
        "n_manifest_c",
        "n_manifest_c_positions",
        "cpg_density_per_100bp",
        "region_rank_by_density",
    ],
    "biomarker_region_cpg": [
        "region_id",
        "site_id",
        "gene_symbol",
        "site_gene_symbols_all",
        "chr",
        "start_pos",
        "cpg_order",
    ],
    "biomarker_region_sequence_score": [
        "region_id",
        "gene_symbol",
        "chr",
        "sequence_region",
        "sequence_start",
        "sequence_end",
        "sequence_available",
        "sequence_error",
        "sequence_length",
        "n_c_sequence",
        "n_g_sequence",
        "n_cg_sequence",
        "n_gcgc",
        "gc_fraction",
        "cg_density_per_100bp",
        "gcgc_density_per_100bp",
        "manifest_cpg_density_per_100bp",
        "sequence_score",
    ],
}


# ============================================================
# DDL — DROP TABLES
# ============================================================

DDL_DROP = """
DROP TABLE IF EXISTS biomarker_region_sequence_score CASCADE;
DROP TABLE IF EXISTS biomarker_region_cpg CASCADE;
DROP TABLE IF EXISTS biomarker_region CASCADE;
DROP TABLE IF EXISTS biomarker_cpg_score CASCADE;
DROP TABLE IF EXISTS cpg_annotation CASCADE;
DROP TABLE IF EXISTS cpg_gene_map CASCADE;
DROP TABLE IF EXISTS expression_correlation CASCADE;
DROP TABLE IF EXISTS cpg_features CASCADE;
DROP TABLE IF EXISTS sample_metadata CASCADE;
DROP TABLE IF EXISTS gene_annotation CASCADE;
DROP TABLE IF EXISTS tumor_types CASCADE;
DROP TABLE IF EXISTS tumor_summary CASCADE;
"""


# ============================================================
# DDL — CREATE TABLES
# ============================================================

TABLES = [
    """
    CREATE TABLE cpg_annotation (
        site_id              TEXT PRIMARY KEY,
        gene                 TEXT,
        chr                  TEXT,
        start_pos            INTEGER,
        end_pos              INTEGER,
        distance_tss         TEXT,
        cgi                  TEXT,
        cg_island            TEXT,
        relation_to_island   TEXT,
        ucsc_refgene_group   TEXT
    );
    """,
    """
    CREATE TABLE cpg_gene_map (
        site_id     TEXT NOT NULL,
        ensembl_id  TEXT,
        gene_symbol TEXT NOT NULL
    );
    """,
    """
    CREATE TABLE expression_correlation (
        site_id      TEXT NOT NULL,
        tumor_type   TEXT NOT NULL,
        sample_class TEXT,
        ensembl_id   TEXT,
        gene_symbol  TEXT,
        spearman_r   REAL,
        pvalue       REAL,
        n_samples    INTEGER
    );
    """,
    """
    CREATE TABLE cpg_features (
        site_id           TEXT PRIMARY KEY,
        leukocyte_median  REAL,
        leukocyte_std     REAL,
        n_samples         INTEGER
    );
    """,
    """
    CREATE TABLE sample_metadata (
        sample_id      TEXT PRIMARY KEY,
        patient_id     TEXT,
        tumor_type     TEXT,
        sample_type    TEXT,
        sample_type_id INTEGER,
        sex            TEXT,
        age            INTEGER,
        tissue_type    TEXT,
        sample_class   TEXT,
        platform       TEXT,
        batch          TEXT
    );
    """,
    """
    CREATE TABLE gene_annotation (
        ensembl_id  TEXT PRIMARY KEY,
        gene_symbol TEXT,
        chr         TEXT,
        start_pos   INTEGER,
        end_pos     INTEGER,
        strand      TEXT,
        biotype     TEXT,
        tss         INTEGER
    );
    """,
    """
    CREATE TABLE tumor_types (
        tumor_type TEXT PRIMARY KEY,
        full_name  TEXT,
        tissue     TEXT
    );
    """,
    """
    CREATE TABLE tumor_summary (
        site_id            TEXT NOT NULL,
        tumor_type         TEXT NOT NULL,
        tumor_median       REAL,
        tumor_std          REAL,
        tumor_n            INTEGER,
        n_tumor            INTEGER,
        normal_median      REAL,
        normal_std         REAL,
        normal_n           INTEGER,
        n_normal           INTEGER,
        pan_tumor_median   REAL,
        pan_tumor_std      REAL,
        pan_tumor_n        INTEGER,
        pan_normal_median  REAL,
        pan_normal_std     REAL,
        pan_normal_n       INTEGER,
        delta_median       REAL,
        hi_index           REAL,
        dispersion_index   REAL,
        PRIMARY KEY (site_id, tumor_type)
    );
    """,
    """
    CREATE TABLE biomarker_cpg_score (
        site_id                       TEXT NOT NULL,
        tumor_type                    TEXT NOT NULL,
        gene_symbol                   TEXT NOT NULL,
        biological_score              REAL,
        delta_score                   REAL,
        normal_low_score              REAL,
        leukocyte_low_score           REAL,
        pancancer_specificity_score   REAL,
        hi_score                      REAL,
        expression_score              REAL,
        sample_support_score          REAL,
        annotation_score              REAL,
        passes_loose_seed             BOOLEAN,
        passes_default_filter         BOOLEAN,
        passes_strict_filter          BOOLEAN,
        PRIMARY KEY (site_id, tumor_type, gene_symbol)
    );
    """,
    """
    CREATE TABLE biomarker_region (
        region_id                 TEXT PRIMARY KEY,
        gene_symbol               TEXT NOT NULL,
        gene_symbols_all          TEXT,
        n_associated_genes        INTEGER,
        chr                       TEXT NOT NULL,
        core_start                INTEGER NOT NULL,
        core_end                  INTEGER NOT NULL,
        core_length               INTEGER,
        browser_start             INTEGER,
        browser_end               INTEGER,
        browser_length            INTEGER,
        flank_bp                  INTEGER,
        n_manifest_cpgs           INTEGER,
        n_manifest_c              INTEGER,
        n_manifest_c_positions    INTEGER,
        cpg_density_per_100bp     REAL,
        region_rank_by_density    INTEGER
    );
    """,
    """
    CREATE TABLE biomarker_region_cpg (
        region_id              TEXT NOT NULL,
        site_id                TEXT NOT NULL,
        gene_symbol            TEXT NOT NULL,
        site_gene_symbols_all  TEXT,
        chr                    TEXT NOT NULL,
        start_pos              INTEGER NOT NULL,
        cpg_order              INTEGER,
        PRIMARY KEY (region_id, site_id, gene_symbol)
    );
    """,
    """
    CREATE TABLE biomarker_region_sequence_score (
        region_id                       TEXT PRIMARY KEY,
        gene_symbol                     TEXT,
        chr                             TEXT,
        sequence_region                 TEXT,
        sequence_start                  INTEGER,
        sequence_end                    INTEGER,
        sequence_available              BOOLEAN,
        sequence_error                  TEXT,
        sequence_length                 INTEGER,
        n_c_sequence                    INTEGER,
        n_g_sequence                    INTEGER,
        n_cg_sequence                   INTEGER,
        n_gcgc                          INTEGER,
        gc_fraction                     REAL,
        cg_density_per_100bp            REAL,
        gcgc_density_per_100bp          REAL,
        manifest_cpg_density_per_100bp  REAL,
        sequence_score                  REAL
    );
    """,
]


# ============================================================
# DDL — INDEXES
# ============================================================

INDEXES = """
CREATE INDEX IF NOT EXISTS idx_tumor_summary_tumor
    ON tumor_summary(tumor_type);
CREATE INDEX IF NOT EXISTS idx_tumor_summary_delta
    ON tumor_summary(delta_median DESC);
CREATE INDEX IF NOT EXISTS idx_tumor_summary_hi
    ON tumor_summary(hi_index DESC);
CREATE INDEX IF NOT EXISTS idx_tumor_summary_dispersion
    ON tumor_summary(dispersion_index DESC);
CREATE INDEX IF NOT EXISTS idx_tumor_summary_site
    ON tumor_summary(site_id);
CREATE INDEX IF NOT EXISTS idx_tumor_summary_site_tumor
    ON tumor_summary(site_id, tumor_type);

CREATE INDEX IF NOT EXISTS idx_cpg_gene_symbol
    ON cpg_gene_map(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_cpg_ensembl
    ON cpg_gene_map(ensembl_id);
CREATE INDEX IF NOT EXISTS idx_cpg_gene_site
    ON cpg_gene_map(site_id);
CREATE INDEX IF NOT EXISTS idx_cpg_gene_site_symbol
    ON cpg_gene_map(site_id, gene_symbol);

CREATE INDEX IF NOT EXISTS idx_expr_corr_site
    ON expression_correlation(site_id);
CREATE INDEX IF NOT EXISTS idx_expr_corr_gene
    ON expression_correlation(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_expr_corr_rho
    ON expression_correlation(spearman_r);
CREATE INDEX IF NOT EXISTS idx_expr_corr_tumor
    ON expression_correlation(tumor_type);
CREATE INDEX IF NOT EXISTS idx_expr_corr_site_tumor_gene
    ON expression_correlation(site_id, tumor_type, gene_symbol);

CREATE INDEX IF NOT EXISTS idx_gene_annotation_symbol
    ON gene_annotation(gene_symbol);

CREATE INDEX IF NOT EXISTS idx_cpg_annotation_chr_start_pos
    ON cpg_annotation(chr, start_pos);
CREATE INDEX IF NOT EXISTS idx_cpg_annotation_gene
    ON cpg_annotation(gene);

CREATE INDEX IF NOT EXISTS idx_cpg_features_site
    ON cpg_features(site_id);

CREATE INDEX IF NOT EXISTS idx_biomarker_cpg_score_tumor_gene
    ON biomarker_cpg_score(tumor_type, gene_symbol);
CREATE INDEX IF NOT EXISTS idx_biomarker_cpg_score_site_tumor_gene
    ON biomarker_cpg_score(site_id, tumor_type, gene_symbol);
CREATE INDEX IF NOT EXISTS idx_biomarker_cpg_score_bio
    ON biomarker_cpg_score(tumor_type, biological_score DESC);

CREATE INDEX IF NOT EXISTS idx_biomarker_region_gene
    ON biomarker_region(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_biomarker_region_genes_all
    ON biomarker_region(gene_symbols_all);
CREATE INDEX IF NOT EXISTS idx_biomarker_region_chr_pos
    ON biomarker_region(chr, core_start, core_end);

CREATE INDEX IF NOT EXISTS idx_biomarker_region_cpg_region
    ON biomarker_region_cpg(region_id);
CREATE INDEX IF NOT EXISTS idx_biomarker_region_cpg_site_gene
    ON biomarker_region_cpg(site_id, gene_symbol);

CREATE INDEX IF NOT EXISTS idx_biomarker_sequence_region
    ON biomarker_region_sequence_score(region_id);
"""


# ============================================================
# HELPERS
# ============================================================

def copy_df_to_table(df: pd.DataFrame, table_name: str, raw_conn) -> None:
    """Insert a DataFrame into PostgreSQL using COPY FROM STDIN."""
    buf = StringIO()
    df.to_csv(buf, index=False, header=False, na_rep="\\N")
    buf.seek(0)

    cols = ", ".join(f'"{c}"' for c in df.columns)

    with raw_conn.cursor() as cur:
        cur.copy_expert(
            f"COPY {table_name} ({cols}) FROM STDIN WITH (FORMAT CSV, NULL '\\N')",
            buf,
        )


def resolve_path(value: Path | Sequence[Path]) -> Path | None:
    """Return the first existing path from a single Path or list of fallback Paths."""
    if isinstance(value, Path):
        return value if value.exists() else None
    for path in value:
        if path.exists():
            return path
    return None


def read_parquet(path: Path) -> pd.DataFrame:
    """Read parquet and print a small summary."""
    df = pd.read_parquet(path)
    print(f"   Rows    : {len(df):,}")
    print(f"   Columns : {list(df.columns)}")
    return df


def safe_numeric(series: pd.Series) -> pd.Series:
    """Convert to numeric while preserving missing values."""
    return pd.to_numeric(series, errors="coerce")


def normalize_gene_symbol(value) -> str:
    """Normalize a gene symbol for CpG-gene and region joins."""
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def split_gene_symbols(value) -> list[str]:
    """Split semicolon/comma-separated gene annotations and preserve order."""
    if pd.isna(value):
        return []
    genes = re.split(r"[;,]", str(value))
    out = []
    for gene in genes:
        gene = normalize_gene_symbol(gene)
        if gene and gene.lower() not in {"nan", "none", "null"}:
            out.append(gene)
    return list(dict.fromkeys(out))


def first_existing_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first existing column among candidates."""
    for col in candidates:
        if col in df.columns:
            return col
    return None


def fill_from_alias(df: pd.DataFrame, target: str, aliases: list[str]) -> pd.DataFrame:
    """Create/fill target from the first available alias column."""
    if target not in df.columns:
        df[target] = pd.NA
    for alias in aliases:
        if alias in df.columns:
            df[target] = df[target].fillna(df[alias])
    return df


def normalize_chrom(series: pd.Series) -> pd.Series:
    """Normalize chromosome labels by removing chr/CHR prefixes."""
    return series.astype(str).str.replace("chr", "", case=False, regex=False).str.strip()


def split_region_universe(universe: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Convert older biomarker_region_universe.parquet into normalized region and region-CpG tables.
    """
    region = universe.copy()
    if "n_manifest_c" not in region.columns and "n_manifest_c_positions" in region.columns:
        region["n_manifest_c"] = region["n_manifest_c_positions"]
    if "n_manifest_c_positions" not in region.columns and "n_manifest_c" in region.columns:
        region["n_manifest_c_positions"] = region["n_manifest_c"]
    if "gene_symbols_all" not in region.columns and "gene_symbol" in region.columns:
        region["gene_symbols_all"] = region["gene_symbol"]
    if "n_associated_genes" not in region.columns and "gene_symbols_all" in region.columns:
        region["n_associated_genes"] = region["gene_symbols_all"].apply(lambda x: len(split_gene_symbols(x)) or pd.NA)

    region_expected = TABLE_COLUMNS["biomarker_region"]
    for col in region_expected:
        if col not in region.columns:
            region[col] = pd.NA

    region = region[region_expected].copy()
    region["gene_symbol"] = region["gene_symbol"].map(normalize_gene_symbol)
    region["chr"] = normalize_chrom(region["chr"])

    for col in [c for c in region_expected if c not in {"region_id", "gene_symbol", "gene_symbols_all", "chr"}]:
        region[col] = safe_numeric(region[col])

    region = region.dropna(subset=["region_id", "gene_symbol", "chr", "core_start", "core_end"])
    region = region.drop_duplicates("region_id")

    region_cpg_columns = TABLE_COLUMNS["biomarker_region_cpg"]
    rows = []
    required = {"region_id", "gene_symbol", "chr", "manifest_site_ids", "manifest_positions"}

    if required.issubset(universe.columns):
        for _, row in universe.iterrows():
            sites = [x for x in str(row.get("manifest_site_ids", "")).split(";") if x]
            positions = [x for x in str(row.get("manifest_positions", "")).split(";") if x]
            for idx, site_id in enumerate(sites):
                pos = positions[idx] if idx < len(positions) else pd.NA
                rows.append(
                    {
                        "region_id": row["region_id"],
                        "site_id": str(site_id),
                        "gene_symbol": normalize_gene_symbol(row["gene_symbol"]),
                        "site_gene_symbols_all": normalize_gene_symbol(row.get("gene_symbols_all", row["gene_symbol"])),
                        "chr": str(row["chr"]).replace("chr", "").replace("CHR", "").strip(),
                        "start_pos": pd.to_numeric(pos, errors="coerce"),
                        "cpg_order": idx + 1,
                    }
                )

    region_cpg = pd.DataFrame(rows, columns=region_cpg_columns)
    if not region_cpg.empty:
        region_cpg["start_pos"] = safe_numeric(region_cpg["start_pos"])
        region_cpg["cpg_order"] = safe_numeric(region_cpg["cpg_order"])
        region_cpg = region_cpg.dropna(subset=["region_id", "site_id", "gene_symbol", "chr", "start_pos"])
        region_cpg = region_cpg.drop_duplicates(["region_id", "site_id", "gene_symbol"])

    return region, region_cpg


def prepare_special_table(df: pd.DataFrame, table_name: str, source_path: Path | None = None) -> pd.DataFrame:
    """Apply table-specific preparation before generic column alignment."""
    out = df.copy()

    if table_name == "cpg_annotation":
        rename = {}
        site_col = first_existing_column(out, ["site_id", "IlmnID", "illumina_id", "probe_id", "Name", "cpg_id"])
        gene_col = first_existing_column(out, ["gene", "gene_symbol", "UCSC_RefGene_Name", "Gene", "gene_name"])
        chr_col = first_existing_column(out, ["chr", "chrom", "chromosome", "CHR", "Chromosome"])
        pos_col = first_existing_column(out, ["start_pos", "position", "pos", "MAPINFO", "start", "Start"])
        end_col = first_existing_column(out, ["end_pos", "end", "End"])
        if site_col and site_col != "site_id":
            rename[site_col] = "site_id"
        if gene_col and gene_col != "gene":
            rename[gene_col] = "gene"
        if chr_col and chr_col != "chr":
            rename[chr_col] = "chr"
        if pos_col and pos_col != "start_pos":
            rename[pos_col] = "start_pos"
        if end_col and end_col != "end_pos":
            rename[end_col] = "end_pos"
        out = out.rename(columns=rename)
        if "chr" in out.columns:
            out["chr"] = normalize_chrom(out["chr"])
        if "end_pos" not in out.columns and "start_pos" in out.columns:
            out["end_pos"] = out["start_pos"]
        if "UCSC_RefGene_Group" in out.columns and "ucsc_refgene_group" not in out.columns:
            out = out.rename(columns={"UCSC_RefGene_Group": "ucsc_refgene_group"})
        if "Relation_to_UCSC_CpG_Island" in out.columns and "relation_to_island" not in out.columns:
            out = out.rename(columns={"Relation_to_UCSC_CpG_Island": "relation_to_island"})

    elif table_name == "cpg_gene_map":
        site_col = first_existing_column(out, ["site_id", "IlmnID", "probe_id", "Name", "cpg_id"])
        gene_col = first_existing_column(out, ["gene_symbol", "gene", "UCSC_RefGene_Name", "Gene", "symbol"])
        ens_col = first_existing_column(out, ["ensembl_id", "ensembl_gene_id", "Ensembl_ID"])
        if site_col and site_col != "site_id":
            out = out.rename(columns={site_col: "site_id"})
        if ens_col and ens_col != "ensembl_id":
            out = out.rename(columns={ens_col: "ensembl_id"})
        if gene_col is None:
            raise ValueError("cpg_gene_map needs a gene_symbol-equivalent column")
        out["gene_symbol_original"] = out[gene_col]
        out["gene_symbol"] = out["gene_symbol_original"].apply(split_gene_symbols)
        out = out.explode("gene_symbol").drop(columns=["gene_symbol_original"])
        out["gene_symbol"] = out["gene_symbol"].map(normalize_gene_symbol)
        out = out[out["gene_symbol"].astype(str).str.len() > 0]

    elif table_name == "expression_correlation":
        rename = {
            "rho": "spearman_r",
            "p_value": "pvalue",
            "pval": "pvalue",
            "n": "n_samples",
        }
        out = out.rename(columns={k: v for k, v in rename.items() if k in out.columns})
        if "gene_symbol" in out.columns:
            out["gene_symbol"] = out["gene_symbol"].fillna("").map(normalize_gene_symbol)
        if "sample_class" not in out.columns:
            out["sample_class"] = "tumor"

    elif table_name == "tumor_summary":
        alias_map = {
            "tumor_n": ["tumor_n", "n_tumor"],
            "n_tumor": ["n_tumor", "tumor_n"],
            "normal_n": ["normal_n", "n_normal"],
            "n_normal": ["n_normal", "normal_n"],
            "pan_tumor_median": ["pan_tumor_median", "panTumor_median", "pantumor_median"],
            "pan_tumor_std": ["pan_tumor_std", "panTumor_std", "pantumor_std"],
            "pan_tumor_n": ["pan_tumor_n", "panTumor_n", "pantumor_n"],
            "pan_normal_median": ["pan_normal_median", "panNormal_median", "pannormal_median"],
            "pan_normal_std": ["pan_normal_std", "panNormal_std", "pannormal_std"],
            "pan_normal_n": ["pan_normal_n", "panNormal_n", "pannormal_n"],
            "hi_index": ["hi_index", "HI_index", "dispersion_index"],
            "dispersion_index": ["dispersion_index", "HI_index", "hi_index"],
        }
        for target, aliases in alias_map.items():
            out = fill_from_alias(out, target, aliases)

    elif table_name == "biomarker_cpg_score":
        rename = {
            "HI_score": "hi_score",
            "passes_default_cpg_filters": "passes_default_filter",
            "passes_strict_cpg_filters": "passes_strict_filter",
            "passes_loose_region_seed": "passes_loose_seed",
        }
        out = out.rename(columns={k: v for k, v in rename.items() if k in out.columns})
        if "gene_symbol" in out.columns:
            out["gene_symbol"] = out["gene_symbol"].map(normalize_gene_symbol)

    elif table_name == "biomarker_region":
        # Current pipeline 14 outputs n_manifest_c_positions; the pages also use n_manifest_c.
        if "n_manifest_c" not in out.columns and "n_manifest_c_positions" in out.columns:
            out["n_manifest_c"] = out["n_manifest_c_positions"]
        if "n_manifest_c_positions" not in out.columns and "n_manifest_c" in out.columns:
            out["n_manifest_c_positions"] = out["n_manifest_c"]
        if "gene_symbol" in out.columns:
            out["gene_symbol"] = out["gene_symbol"].map(normalize_gene_symbol)
        if "gene_symbols_all" not in out.columns and "gene_symbol" in out.columns:
            out["gene_symbols_all"] = out["gene_symbol"]
        if "gene_symbols_all" in out.columns:
            out["gene_symbols_all"] = out["gene_symbols_all"].apply(lambda x: ";".join(split_gene_symbols(x)))
        if "n_associated_genes" not in out.columns and "gene_symbols_all" in out.columns:
            out["n_associated_genes"] = out["gene_symbols_all"].apply(lambda x: len(split_gene_symbols(x)) or pd.NA)
        if "chr" in out.columns:
            out["chr"] = normalize_chrom(out["chr"])

    elif table_name == "biomarker_region_cpg":
        if "gene_symbol" in out.columns:
            out["gene_symbol"] = out["gene_symbol"].map(normalize_gene_symbol)
        if "site_gene_symbols_all" not in out.columns and "gene_symbol" in out.columns:
            out["site_gene_symbols_all"] = out["gene_symbol"]
        if "site_gene_symbols_all" in out.columns:
            out["site_gene_symbols_all"] = out["site_gene_symbols_all"].apply(lambda x: ";".join(split_gene_symbols(x)))
        if "chr" in out.columns:
            out["chr"] = normalize_chrom(out["chr"])

    elif table_name == "biomarker_region_sequence_score":
        if "gene_symbol" in out.columns:
            out["gene_symbol"] = out["gene_symbol"].map(normalize_gene_symbol)
        if "chr" in out.columns:
            out["chr"] = normalize_chrom(out["chr"])
        if "manifest_cpg_density_per_100bp" not in out.columns and "cpg_density_manifest_per_100bp" in out.columns:
            out["manifest_cpg_density_per_100bp"] = out["cpg_density_manifest_per_100bp"]

    return out


def align_columns(df: pd.DataFrame, table_name: str) -> pd.DataFrame:
    """
    Reorder columns to match the PostgreSQL table.
    Missing optional columns are filled with NA.
    Extra columns are ignored.
    """
    expected = TABLE_COLUMNS[table_name]

    for col in expected:
        if col not in df.columns:
            df[col] = pd.NA
            print(f"   ⚠ Missing column filled with NULL: {col}")

    extra = [col for col in df.columns if col not in expected]
    if extra:
        print(f"   ⚠ Extra columns ignored: {extra[:20]}" + (" ..." if len(extra) > 20 else ""))

    out = df[expected].copy()

    # Type cleanup for commonly numeric columns.
    int_like = {
        "start_pos", "end_pos", "tumor_n", "normal_n", "n_tumor", "n_normal",
        "pan_tumor_n", "pan_normal_n", "n_samples",
        "n_leukocyte", "sample_type_id", "age", "tss", "core_start", "core_end",
        "core_length", "browser_start", "browser_end", "browser_length", "flank_bp",
        "n_manifest_cpgs", "n_manifest_c", "n_manifest_c_positions", "n_associated_genes",
        "region_rank_by_density", "cpg_order", "sequence_start", "sequence_end",
        "sequence_length", "n_c_sequence", "n_g_sequence", "n_cg_sequence", "n_gcgc",
    }
    float_like = {
        "tumor_median", "tumor_std", "normal_median", "normal_std", "pan_tumor_median",
        "pan_tumor_std", "pan_normal_median", "pan_normal_std", "delta_median",
        "hi_index", "dispersion_index", "leukocyte_median", "leukocyte_std", "spearman_r",
        "pvalue", "biological_score", "delta_score", "normal_low_score",
        "leukocyte_low_score", "pancancer_specificity_score", "hi_score", "expression_score",
        "sample_support_score", "annotation_score", "cpg_density_per_100bp", "gc_fraction",
        "cg_density_per_100bp", "gcgc_density_per_100bp", "manifest_cpg_density_per_100bp",
        "sequence_score",
    }
    bool_like = {"passes_loose_seed", "passes_default_filter", "passes_strict_filter", "sequence_available"}

    for col in int_like.intersection(out.columns):
        out[col] = safe_numeric(out[col]).astype("Int64")
    for col in float_like.intersection(out.columns):
        out[col] = safe_numeric(out[col])
    for col in bool_like.intersection(out.columns):
        out[col] = out[col].fillna(False).astype(bool)

    return out


def deduplicate(df: pd.DataFrame, table_name: str) -> pd.DataFrame:
    """Apply table-specific deduplication rules."""
    rules = {
        "cpg_gene_map": ["site_id", "ensembl_id", "gene_symbol"],
        "expression_correlation": ["site_id", "tumor_type", "sample_class", "ensembl_id", "gene_symbol"],
        "tumor_summary": ["site_id", "tumor_type"],
        "sample_metadata": ["sample_id"],
        "gene_annotation": ["ensembl_id"],
        "cpg_annotation": ["site_id"],
        "cpg_features": ["site_id"],
        "tumor_types": ["tumor_type"],
        "biomarker_cpg_score": ["site_id", "tumor_type", "gene_symbol"],
        "biomarker_region": ["region_id"],
        "biomarker_region_cpg": ["region_id", "site_id", "gene_symbol"],
        "biomarker_region_sequence_score": ["region_id"],
    }

    subset = [col for col in rules.get(table_name, []) if col in df.columns]
    if not subset:
        return df

    before = len(df)
    df = df.drop_duplicates(subset=subset)
    removed = before - len(df)

    if removed:
        print(f"   Dedup   : {removed:,} duplicates removed using {subset}")

    return df


def should_load_table(table_name: str, args: argparse.Namespace) -> bool:
    if args.skip_base and table_name in BASE_TABLES:
        return False
    if args.skip_biomarkers and table_name in BIOMARKER_TABLES:
        return False
    return True


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create/load PostgreSQL tables for methylation app plus region biomarker model."
    )
    parser.add_argument("--skip-base", action="store_true", help="Do not drop/create/load base app tables.")
    parser.add_argument("--skip-biomarkers", action="store_true", help="Do not drop/create/load biomarker region tables.")
    parser.add_argument("--schema-only", action="store_true", help="Create schema and indexes only; do not import parquet files.")
    parser.add_argument("--no-drop", action="store_true", help="Do not drop existing tables before creating schema.")
    return parser.parse_args()


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    args = parse_args()
    engine = get_engine()

    print("\n🏗  Creating PostgreSQL schema...\n")

    with engine.begin() as conn:
        if not args.no_drop:
            conn.execute(text(DDL_DROP))
            print("✓ Old project/biomarker tables dropped")
        else:
            print("✓ Keeping existing tables because --no-drop was used")

        for i, sql in enumerate(TABLES, 1):
            try:
                conn.execute(text(sql))
                print(f"✓ Table {i}/{len(TABLES)} created")
            except Exception as e:
                print("\n================ ERROR ================")
                print(f"  Table #{i}")
                print(f"  Type : {type(e).__name__}")
                print(f"  Error: {e.orig if hasattr(e, 'orig') else e}")
                print("=======================================\n")
                traceback.print_exc(limit=2)
                raise SystemExit(1)

        conn.execute(text(INDEXES))
        print("✓ Indexes created")

    if args.schema_only:
        engine.dispose()
        print("\n✅ Schema ready. No data imported because --schema-only was used.")
        return

    print("\n🚀 Importing parquet files...\n")

    raw_conn = engine.raw_connection()
    imported = []
    skipped = []
    failed = []

    try:
        region_universe_cache: pd.DataFrame | None = None

        for table_name, parquet_candidate in FILES.items():
            if not should_load_table(table_name, args):
                print(f"⚠  Skipping {table_name} — disabled by CLI flag")
                skipped.append(table_name)
                continue

            parquet_path = resolve_path(parquet_candidate)
            if parquet_path is None:
                print(f"⚠  Skipping {table_name} — file not found: {parquet_candidate}")
                skipped.append(table_name)
                continue

            print(f"\n📦 {table_name}")
            print(f"   File : {parquet_path}")

            try:
                # If using the older region universe output, split it into both normalized tables.
                if table_name in {"biomarker_region", "biomarker_region_cpg"} and parquet_path.name == "biomarker_region_universe.parquet":
                    if region_universe_cache is None:
                        region_universe_cache = read_parquet(parquet_path)
                    region_df, region_cpg_df = split_region_universe(region_universe_cache)
                    df = region_df if table_name == "biomarker_region" else region_cpg_df
                    print(f"   Converted from region universe: shape={df.shape}")
                else:
                    df = read_parquet(parquet_path)
                    df = prepare_special_table(df, table_name, source_path=parquet_path)

                df = align_columns(df, table_name)
                df = deduplicate(df, table_name)

                if df.empty:
                    print(f"   ⚠ Empty after preparation. Skipping import.")
                    skipped.append(table_name)
                    continue

                try:
                    copy_df_to_table(df, table_name, raw_conn)
                    raw_conn.commit()
                    imported.append(table_name)
                    print(f"   ✓ Imported {len(df):,} rows")
                except Exception as e:
                    raw_conn.rollback()
                    failed.append(table_name)
                    print(f"\n   ✗ ERROR importing {table_name}")
                    err_lines = str(e).splitlines()
                    for line in err_lines[:12]:
                        print(f"   {line}")
                    if len(err_lines) > 12:
                        print(f"   ... ({len(err_lines) - 12} lines hidden)")
                    print()

            except Exception as e:
                failed.append(table_name)
                print(f"\n   ✗ ERROR preparing {table_name}: {type(e).__name__}: {e}")
                traceback.print_exc(limit=2)

    finally:
        raw_conn.close()

    print("\n📊 Running ANALYZE...")
    with engine.begin() as conn:
        conn.execute(text("ANALYZE"))
    print("✓ ANALYZE completed")

    engine.dispose()

    print("\n================ SUMMARY ================")
    print(f"Imported: {imported}")
    print(f"Skipped : {skipped}")
    print(f"Failed  : {failed}")
    print("=========================================\n")

    if failed:
        raise SystemExit("Some tables failed to import. Check messages above.")

    print("✅ PostgreSQL database ready.")


if __name__ == "__main__":
    main()
