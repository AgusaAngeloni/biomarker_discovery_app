from io import StringIO
from pathlib import Path

import pandas as pd
from sqlalchemy import create_engine

ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = ROOT / "data" / "raw"

POSTGRES_USER = "neondb_owner"
POSTGRES_PASSWORD = "methy_2026_EpiLiquiD"

POSTGRES_HOST = "ep-super-flower-aq09lnkw-pooler.c-8.us-east-1.aws.neon.tech"
POSTGRES_PORT = 5432
POSTGRES_DB = "neondb"

DB_URL = (
    f"postgresql+psycopg2://"
    f"{POSTGRES_USER}:{POSTGRES_PASSWORD}"
    f"@{POSTGRES_HOST}:{POSTGRES_PORT}/{POSTGRES_DB}"
    "?sslmode=require"
)

FILES = [
    ("expression_correlation", DATA_DIR / "expression_correlation.parquet")
    
]

def copy_df_to_table(df, table_name, raw_conn):

    buf = StringIO()

    df.to_csv(
        buf,
        index=False,
        header=False,
        na_rep="\\N"
    )

    buf.seek(0)

    cols = ", ".join(f'"{c}"' for c in df.columns)

    with raw_conn.cursor() as cur:
        cur.copy_expert(
            f"""
            COPY {table_name} ({cols})
            FROM STDIN
            WITH (FORMAT CSV, NULL '\\N')
            """,
            buf,
        )

print("\n🔌 Conectando a Neon...\n")

engine = create_engine(DB_URL)
raw_conn = engine.raw_connection()

try:

    for table_name, parquet_file in FILES:

        print(f"\n📦 {table_name}")

        df = pd.read_parquet(parquet_file)

        with raw_conn.cursor() as cur:
            cur.execute(f"TRUNCATE TABLE {table_name}")

        raw_conn.commit()

        copy_df_to_table(
            df,
            table_name,
            raw_conn
        )

        raw_conn.commit()

        print(f"✓ {table_name}: {len(df):,} filas")

finally:

    raw_conn.close()
    engine.dispose()

print("\n✅ Carga terminada")
