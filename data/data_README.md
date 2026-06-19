# Data

This folder stores raw and derived data files used by the preprocessing pipeline and PostgreSQL import step.

## Expected structure

```text
data/
├── raw/
├── summary/
├── correlations_tumor/
└── geo/
    └── IDATS/
```

## Main derived files expected by `build_postgres.py`

```text
data/raw/manifest_clean.parquet
data/raw/cpg_gene_map.parquet
data/raw/expression_correlation.parquet
data/raw/cpg_features.parquet
data/raw/pheno_clean.parquet
data/raw/gene_annotation.parquet
data/raw/tumor_types.parquet
data/raw/methylation_summary.parquet
```

## GitHub policy

Large raw and derived files should not be committed to GitHub.

Recommended `.gitignore` entries:

```text
data/raw/*.parquet
data/raw/*.tsv
data/raw/*.csv
data/summary/
data/correlations_tumor/
data/geo/
```

For reproducibility, document the source of each raw file and provide scripts to regenerate derived outputs.
