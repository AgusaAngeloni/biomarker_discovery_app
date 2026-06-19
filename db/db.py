from __future__ import annotations

import os
from sqlalchemy import create_engine


def get_database_url() -> str:
    """
    Get database URL from Streamlit secrets or environment variable.

    Priority:
    1. st.secrets["database"]["url"]
    2. DATABASE_URL environment variable
    """

    try:
        import streamlit as st

        if "database" in st.secrets and "url" in st.secrets["database"]:
            return st.secrets["database"]["url"]

    except Exception:
        pass

    database_url = os.getenv("DATABASE_URL")

    if database_url:
        return database_url

    raise RuntimeError(
        "Database URL not found. Define it in "
        ".streamlit/secrets.toml as [database] url = '...', "
        "or export DATABASE_URL."
    )


def get_engine():
    database_url = get_database_url()

    return create_engine(
        database_url,
        pool_pre_ping=True,
    )
