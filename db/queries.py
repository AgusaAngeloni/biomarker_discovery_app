#queries.py

from sqlalchemy import create_engine, text
import pandas as pd


DB_URL = (
    "postgresql+psycopg2://"
    "postgres:methy20@localhost:5432/db_methylation"
)

engine = create_engine(DB_URL)

def run_query(query, params=None):

    with engine.connect() as conn:

        return pd.read_sql(
            text(query),
            conn,
            params=params
        )
