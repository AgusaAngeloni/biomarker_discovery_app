# Database

This folder contains PostgreSQL utilities and SQL-related files used by the Streamlit app.

## Local PostgreSQL setup

The database can be created manually:

```bash
createdb db_methylation
```

or automatically by:

```bash
python pipelines/build_postgres.py
```

## Environment variables for local database build

```bash
export POSTGRES_USER=postgres
export POSTGRES_PASSWORD=your_password
export POSTGRES_HOST=localhost
export POSTGRES_PORT=5432
export POSTGRES_DB=db_methylation
```

Then run:

```bash
python pipelines/build_postgres.py
```

## Test local connection

```bash
psql -h localhost -U postgres -d db_methylation
```

Useful PostgreSQL commands:

```sql
\dt
SELECT COUNT(*) FROM tumor_summary;
SELECT COUNT(*) FROM cpg_annotation;
SELECT COUNT(*) FROM expression_correlation;
```

## Streamlit connection

Streamlit should read the database URL from `.streamlit/secrets.toml`:

```toml
[database]
url = "postgresql+psycopg2://postgres:your_password@localhost:5432/db_methylation"
```

A typical `db/database.py` file can expose the SQLAlchemy engine:

```python
import streamlit as st
from sqlalchemy import create_engine

@st.cache_resource
def get_engine():
    return create_engine(st.secrets["database"]["url"])
```

## Main tables

### `cpg_annotation`

CpG probe annotation, including genomic position and CpG island context.

### `cpg_gene_map`

Mapping between CpG probes and genes.

### `expression_correlation`

Methylation-expression correlation per CpG and tumor type.

### `cpg_features`

CpG-level features such as leukocyte methylation background.

### `sample_metadata`

Cleaned sample metadata.

### `gene_annotation`

Gene coordinates, strand, biotype, and TSS information.

### `tumor_types`

Tumor type dictionary.

### `tumor_summary`

Tumor vs normal methylation summary by CpG and tumor type.
