# Database README

This folder contains the lightweight database access layer used by the Streamlit app.

The PostgreSQL schema is created and populated by:

```text
pipelines/build_postgres.py
```

The files in this folder are used at app runtime to connect to the database and execute SQL queries.

---

## Files

```text
db/
├── db.py
├── queries.py
└── db_README.md
```

---

## `db.py`

Defines the database connection logic.

The app checks the database URL in this order:

1. Streamlit secrets:

```toml
[database]
url = "postgresql+psycopg2://user:password@host:5432/dbname"
```

2. Environment variable:

```bash
export DATABASE_URL="postgresql+psycopg2://user:password@host:5432/dbname"
```

If neither is available, the app raises an error asking for a database URL.

---

## `queries.py`

Defines the helper used by Streamlit pages:

```python
run_query(query: str, params: dict | None = None) -> pandas.DataFrame
```

The helper:

- gets the SQLAlchemy engine from `db.py`;
- executes parameterized SQL;
- returns the result as a pandas DataFrame.

Pages should use parameterized queries instead of string interpolation for user-selected values.

---

## Expected PostgreSQL tables

### Base app tables

```text
cpg_annotation
cpg_gene_map
expression_correlation
cpg_features
sample_metadata
gene_annotation
tumor_types
tumor_summary
```

### Region / biomarker tables

```text
biomarker_cpg_score
biomarker_region
biomarker_region_cpg
biomarker_region_sequence_score
```

---

## Important column conventions

The preferred standardized names are:

```text
tumor_summary.hi_index
tumor_summary.pan_tumor_median
tumor_summary.pan_normal_median
```

The app may include compatibility logic for older aliases such as:

```text
dispersion_index
HI_index
pantumor_median
pannormal_median
panTumor_median
panNormal_median
```

However, new database builds should use the standardized names.

---

## Building the database

Run from the project root.

### Local PostgreSQL

```bash
export POSTGRES_USER=postgres
export POSTGRES_PASSWORD=your_password
export POSTGRES_HOST=localhost
export POSTGRES_PORT=5432
export POSTGRES_DB=db_methylation

python pipelines/build_postgres.py
```

### Remote PostgreSQL

```bash
export DATABASE_URL="postgresql+psycopg2://user:password@host:5432/dbname"

python pipelines/build_postgres.py
```

When `DATABASE_URL` is provided, the loader connects directly to that database and does not attempt to create a local database.

---

## Runtime configuration for Streamlit

For local development, either export `DATABASE_URL`:

```bash
export DATABASE_URL="postgresql+psycopg2://user:password@localhost:5432/db_methylation"
streamlit run Main.py
```

or create `.streamlit/secrets.toml`:

```toml
[database]
url = "postgresql+psycopg2://user:password@localhost:5432/db_methylation"
```

Do not commit `.streamlit/secrets.toml`.

---

## Query design notes

- PostgreSQL should be treated as the source of truth for the app.
- Pages should filter in SQL before loading results into pandas.
- Region pages should preserve `biomarker_region.region_id`.
- CpG-to-gene joins should use `cpg_gene_map` rather than parsing compound gene annotations from the manifest.
- Expression joins should preserve the gene context when possible:

```text
site_id + tumor_type + gene_symbol
```

or:

```text
site_id + tumor_type + ensembl_id
```

This avoids assigning a CpG-expression correlation to the wrong gene when a CpG maps to multiple genes.

---

## Schema maintenance

If the database model changes, update together:

```text
pipelines/build_postgres.py
schema.sql
schema_neon.sql
README.md
pages/pages_README.md
pipelines/pipelines_README.md
db/db_README.md
```

The most important rule is that the Streamlit queries, PostgreSQL schema, and Parquet outputs must agree on table and column names.
