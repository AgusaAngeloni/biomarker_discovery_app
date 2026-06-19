"""
build_postgres.py
───────────────────────────────────────────────────────────────
Create the PostgreSQL schema and import derived .parquet tables.

This script is intended for reproducible local/database setup.
It can be used in two ways:

1) Local PostgreSQL using separate environment variables:

    export POSTGRES_USER=postgres
    export POSTGRES_PASSWORD=your_password
    export POSTGRES_HOST=localhost
    export POSTGRES_PORT=5432
    export POSTGRES_DB=db_methylation

    python pipelines/build_postgres.py

2) Direct database URL, useful for remote PostgreSQL providers:

    export DATABASE_URL="postgresql+psycopg2://user:password@host:5432/dbname"
    python pipelines/build_postgres.py

Notes:
    - If DATABASE_URL is provided, the script connects directly to that DB.
      It does not try to create the database itself.
    - If DATABASE_URL is not provided, the script builds the URL from POSTGRES_* variables
      and attempts to create POSTGRES_DB if it does not exist.
    - Existing project tables are dropped and recreated.
    - Large input files should exist under data/raw/.

Requirements:
    pip install pandas pyarrow sqlalchemy psycopg2-binary
"""

from __future__ import annotations

from io import StringIO
from pathlib import Path
import os
import traceback
from typing import Dict, List

import pandas as pd
from sqlalchemy import create_engine, text
from sqlalchemy.engine import URL


# ============================================================
# PATHS
# ============================================================

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data" / "raw"


# ============================================================
# DATABASE CONFIG
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

    If DATABASE_URL is present, connect directly to it.
    Otherwise, create POSTGRES_DB if needed and connect to it.
    """
    if DATABASE_URL:
        print("\n🔗 Using DATABASE_URL")
        return create_engine(DATABASE_URL)

    print("\n🔍 Checking local PostgreSQL database...\n")

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

FILES = {
    "cpg_annotation": DATA_DIR / "manifest_clean.parquet",
    "cpg_gene_map": DATA_DIR / "cpg_gene_map.parquet",
    "expression_correlation": DATA_DIR / "expression_correlation.parquet",
    "cpg_features": DATA_DIR / "cpg_features.parquet",
    "sample_metadata": DATA_DIR / "pheno_clean.parquet",
    "gene_annotation": DATA_DIR / "gene_annotation.parquet",
    "tumor_types": DATA_DIR / "tumor_types.parquet",
    "tumor_summary": DATA_DIR / "methylation_summary.parquet",

}


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
    "tumor_summary": [
        "site_id",
        "tumor_type",
        "tumor_median",
        "tumor_std",
        "tumor_n",
        "normal_median",
        "normal_std",
        "normal_n",
        "pantumor_median",
        "pantumor_std",
        "pantumor_n",
        "pannormal_median",
        "pannormal_std",
        "pannormal_n",
        "delta_median",
        "hi_index",
    ],
}


# ============================================================
# DDL — DROP TABLES
# ============================================================

DDL_DROP = """
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
        site_id       TEXT PRIMARY KEY,
        gene          TEXT,
        chr           TEXT,
        start_pos     INTEGER,
        end_pos       INTEGER,
        distance_tss  TEXT,
        cgi           TEXT,
        cg_island     TEXT
    );
    """,
    """
    CREATE TABLE cpg_gene_map (
        site_id     TEXT NOT NULL,
        ensembl_id  TEXT NOT NULL,
        gene_symbol TEXT
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
        site_id           TEXT NOT NULL,
        tumor_type        TEXT NOT NULL,
        tumor_median      REAL,
        tumor_std         REAL,
        tumor_n           INTEGER,  
        normal_median     REAL,
        normal_std        REAL,
        normal_n          INTEGER,  
        panTumor_median    REAL,
        panTumor_std      REAL,
        panTumor_n        INTEGER,  
        panNormal_median  REAL,
        panNormal_std     REAL,
        panNormal_n       INTEGER,  
        delta_median      REAL,
        hi_index          REAL,
        PRIMARY KEY (site_id, tumor_type)
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
CREATE INDEX IF NOT EXISTS idx_tumor_summary_HI
    ON tumor_summary(hi_index DESC);
CREATE INDEX IF NOT EXISTS idx_tumor_summary_site
    ON tumor_summary(site_id);

CREATE INDEX IF NOT EXISTS idx_cpg_gene_symbol
    ON cpg_gene_map(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_cpg_ensembl
    ON cpg_gene_map(ensembl_id);
CREATE INDEX IF NOT EXISTS idx_cpg_gene_site
    ON cpg_gene_map(site_id);

CREATE INDEX IF NOT EXISTS idx_expr_corr_site
    ON expression_correlation(site_id);
CREATE INDEX IF NOT EXISTS idx_expr_corr_gene
    ON expression_correlation(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_expr_corr_rho
    ON expression_correlation(spearman_r);
CREATE INDEX IF NOT EXISTS idx_expr_corr_tumor
    ON expression_correlation(tumor_type);

CREATE INDEX IF NOT EXISTS idx_gene_annotation_symbol
    ON gene_annotation(gene_symbol);

CREATE INDEX IF NOT EXISTS idx_cpg_annotation_chr_start_pos
    ON cpg_annotation(chr, start_pos);
CREATE INDEX IF NOT EXISTS idx_cpg_annotation_gene
    ON cpg_annotation(gene);
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


def read_parquet(path: Path) -> pd.DataFrame:
    """Read parquet and print a small summary."""
    df = pd.read_parquet(path)
    print(f"   Rows    : {len(df):,}")
    print(f"   Columns : {list(df.columns)}")
    return df


def align_columns(df: pd.DataFrame, table_name: str) -> pd.DataFrame:
    """
    Reorder columns to match the PostgreSQL table.
    Missing optional columns are filled with NA.
    Extra columns are ignored.
    """
    if table_name == "tumor_summary":
        df = df.rename(columns={
            "panTumor_median": "pantumor_median",
            "panTumor_std": "pantumor_std",
            "panTumor_n": "pantumor_n",
            "panNormal_median": "pannormal_median",
            "panNormal_std": "pannormal_std",
            "panNormal_n": "pannormal_n",
            "HI_index": "hi_index",
        })
    
    expected = TABLE_COLUMNS[table_name]

    for col in expected:
        if col not in df.columns:
            df[col] = pd.NA
            print(f"   ⚠ Missing column filled with NULL: {col}")

    extra = [col for col in df.columns if col not in expected]
    if extra:
        print(f"   ⚠ Extra columns ignored: {extra}")

    return df[expected]


def deduplicate(df: pd.DataFrame, table_name: str) -> pd.DataFrame:
    """Apply table-specific deduplication rules."""
    rules = {
        "cpg_gene_map": ["site_id", "ensembl_id"],
        "expression_correlation": ["site_id", "tumor_type", "sample_class", "ensembl_id"],
        "tumor_summary": ["site_id", "tumor_type"],
        "sample_metadata": ["sample_id"],
        "gene_annotation": ["ensembl_id"],
        "cpg_annotation": ["site_id"],
        "cpg_features": ["site_id"],
        "tumor_types": ["tumor_type"],
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


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    engine = get_engine()

    print("\n🏗  Creating PostgreSQL schema...\n")

    with engine.begin() as conn:
        conn.execute(text(DDL_DROP))
        print("✓ Old project tables dropped")

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

    print("\n🚀 Importing parquet files...\n")

    raw_conn = engine.raw_connection()
    imported = []
    skipped = []
    failed = []

    try:
        for table_name, parquet_path in FILES.items():
            if not parquet_path.exists():
                print(f"⚠  Skipping {table_name} — file not found: {parquet_path}")
                skipped.append(table_name)
                continue

            print(f"\n📦 {table_name}")
            print(f"   File : {parquet_path}")

            df = read_parquet(parquet_path)
            df = align_columns(df, table_name)
            df = deduplicate(df, table_name)

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
