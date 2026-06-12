"""
build_postgres.py
───────────────────────────────────────────────────────────────
Crear esquema PostgreSQL e importar tablas .parquet

Requisitos:
    pip install pandas pyarrow sqlalchemy psycopg2-binary

Uso:
    python build_postgres.py

Asegurarse de:
    - PostgreSQL corriendo
    - Variables de conexión correctas
───────────────────────────────────────────────────────────────
"""

from io import StringIO
from pathlib import Path
import traceback

import pandas as pd
from sqlalchemy import create_engine, text

ROOT = Path(__file__).resolve().parent.parent

# ============================================================
# CONFIG
# ============================================================

POSTGRES_USER     = "postgres"
POSTGRES_PASSWORD = "methy20"
POSTGRES_HOST     = "localhost"
POSTGRES_PORT     = 5432
POSTGRES_DB       = "db_methylation"

DATA_DIR = ROOT / "data" / "raw"

# ============================================================
# PARQUETS
# ============================================================

FILES = {
    "cpg_annotation":         DATA_DIR / "manifest_clean.parquet",
    "cpg_gene_map":           DATA_DIR / "cpg_gene_map.parquet",
    "expression_correlation": DATA_DIR / "expression_correlation.parquet",
    "cpg_features":           DATA_DIR / "cpg_features.parquet",
    "sample_metadata":        DATA_DIR / "pheno_clean.parquet",
    "gene_annotation":        DATA_DIR / "gene_annotation.parquet",
    "tumor_types":            DATA_DIR / "tumor_types.parquet",
    "tumor_summary":          DATA_DIR / "methylation_summary.parquet",
}

# ============================================================
# HELPERS
# ============================================================

def copy_df_to_table(df: pd.DataFrame, table_name: str, raw_conn) -> None:
    """
    Inserta un DataFrame en PostgreSQL usando COPY FROM STDIN (CSV).
    Sin límite de parámetros, mucho más rápido que INSERT multi-row.
    """
    buf = StringIO()
    df.to_csv(buf, index=False, header=False, na_rep="\\N")
    buf.seek(0)

    cols = ", ".join(f'"{c}"' for c in df.columns)

    with raw_conn.cursor() as cur:
        cur.copy_expert(
            f"COPY {table_name} ({cols}) FROM STDIN WITH (FORMAT CSV, NULL '\\N')",
            buf,
        )


def peek_parquet(path: Path) -> pd.DataFrame:
    """Lee el parquet y muestra un resumen de columnas y tipos."""
    df = pd.read_parquet(path)
    print(f"   Rows    : {len(df):,}")
    print(f"   Columns : {list(df.columns)}")
    return df


# ============================================================
# CONNECT — check / create database
# ============================================================

admin_url = (
    f"postgresql+psycopg2://"
    f"{POSTGRES_USER}:{POSTGRES_PASSWORD}"
    f"@{POSTGRES_HOST}:{POSTGRES_PORT}/postgres"
)

admin_engine = create_engine(admin_url)

print("\n🔍 Checking database...\n")

with admin_engine.connect() as conn:
    conn = conn.execution_options(isolation_level="AUTOCOMMIT")

    exists = conn.execute(
        text("SELECT 1 FROM pg_database WHERE datname = :dbname"),
        {"dbname": POSTGRES_DB},
    ).scalar()

    if not exists:
        print(f"   Creating: {POSTGRES_DB}")
        conn.execute(text(f'CREATE DATABASE "{POSTGRES_DB}"'))
        print("   ✓ Database created")
    else:
        print("   ✓ Database already exists")

admin_engine.dispose()

# ============================================================
# ENGINE
# ============================================================

DB_URL = (
    f"postgresql+psycopg2://"
    f"{POSTGRES_USER}:{POSTGRES_PASSWORD}"
    f"@{POSTGRES_HOST}:{POSTGRES_PORT}/{POSTGRES_DB}"
)

engine = create_engine(DB_URL)

# ============================================================
# DDL — DROP
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
        start_pos         INTEGER,
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
        -- sample_class TEXT NOT NULL,
        -- ensembl_id   TEXT NOT NULL,
        -- gene_symbol  TEXT, 
        spearman_r   REAL,
        -- pvalue       REAL,
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
        site_id          TEXT NOT NULL,
        tumor_type       TEXT NOT NULL,
        tumor_median     REAL,
        tumor_std        REAL,
        normal_median    REAL,
        normal_std       REAL,
        pan_tumor_median  REAL,
        pan_tumor_std     REAL,
        pan_normal_median REAL,
        pan_normal_std    REAL,
        delta_median     REAL,
        dispersion_index REAL,
        n_tumor          INTEGER,
        n_normal         INTEGER,
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
--CREATE INDEX IF NOT EXISTS idx_tumor_summary_dispersion
  --  ON tumor_summary(dispersion_index DESC);

CREATE INDEX IF NOT EXISTS idx_cpg_gene_symbol
    ON cpg_gene_map(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_cpg_ensembl
    ON cpg_gene_map(ensembl_id);

CREATE INDEX IF NOT EXISTS idx_expr_corr_site
    ON expression_correlation(site_id);
-- CREATE INDEX IF NOT EXISTS idx_expr_corr_gene
   -- ON expression_correlation(gene_symbol);
--CREATE INDEX IF NOT EXISTS idx_expr_corr_rho
  --  ON expression_correlation(spearman_r);
CREATE INDEX IF NOT EXISTS idx_expr_corr_tumor
    ON expression_correlation(tumor_type);

CREATE INDEX IF NOT EXISTS idx_gene_annotation_symbol
    ON gene_annotation(gene_symbol);

CREATE INDEX IF NOT EXISTS idx_cpg_annotation_chr_start_pos
    ON cpg_annotation(chr, start_pos);
"""

# ============================================================
# CREATE SCHEMA
# ============================================================

print("\n🏗  Creating PostgreSQL schema...\n")

with engine.begin() as conn:

    conn.execute(text(DDL_DROP))
    print("✓ Old tables dropped")

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

# ============================================================
# IMPORT PARQUETS  —  via COPY (sin límite de parámetros)
# ============================================================

print("\n🚀 Importing parquet files...\n")

raw_conn = engine.raw_connection()
filed = []
skipped = []

try:
    for table_name, parquet_path in FILES.items():

        if not parquet_path.exists():
            print(f"⚠  Skipping (file not found): {parquet_path}")
            skipped.append(table_name)            
            continue

        print(f"\n📦 {table_name}")
        print(f"   File : {parquet_path}")

        df = peek_parquet(parquet_path)
        
        if table_name == "cpg_gene_map":
            before = len(df)
            df = df.drop_duplicates(subset=["site_id","ensembl_id"])
            print(f"    Dedup   :{before - len(df):,}duplicates removed")
        if table_name == "expression_correlation":
            before = len(df)
            df = df.drop_duplicates(subset=["site_id","tumor_type"])
            print(f"    Dedup   :{before - len(df):,}duplicates removed")

        try:
            copy_df_to_table(df, table_name, raw_conn)
            raw_conn.commit()
            print(f"   ✓ Imported {len(df):,} rows")

        except Exception as e:
            raw_conn.rollback()
            print(f"\n   ✗ ERROR importing {table_name}")
            err_lines = str(e).splitlines()
            for line in err_lines[:10]:
                print(f"   {line}")
            if len(err_lines) > 10:
                print(f"   ... ({len(err_lines) - 10} lines hidden)")
            print()
            continue

finally:
    raw_conn.close()

engine.dispose()
print("\n✅ PostgreSQL database ready.")
