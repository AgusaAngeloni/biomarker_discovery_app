# db/queries.py

from __future__ import annotations

import pandas as pd
from sqlalchemy import text

from db.db import get_engine


def run_query(query: str, params: dict | None = None) -> pd.DataFrame:
    """
    Run a SQL query and return the result as a pandas DataFrame.
    """

    engine = get_engine()

    with engine.connect() as conn:
        return pd.read_sql(
            text(query),
            conn,
            params=params,
        )
