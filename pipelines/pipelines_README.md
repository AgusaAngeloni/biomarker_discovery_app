# Pipelines

This folder contains the preprocessing workflow used to generate the methylation biomarker database.

The scripts convert raw methylation, phenotype, gene annotation, CpG annotation, expression, and leukocyte methylation data into derived `.parquet` files that are later imported into PostgreSQL.

## Recommended execution order

### 1. Input data preparation

```bash
python pipelines/01_clean_manifiesto.py
python pipelines/02_clean_phenotype.py
python pipelines/03_methy_download_clean_v2.py
```

Expected outputs include:

```text
data/raw/manifest_clean.parquet
data/raw/pheno_clean.parquet
```

### 2. Methylation summaries

```bash
python pipelines/04_methy_summary_v2.py
python pipelines/05_merge_methy_summary.py
```

Expected output:

```text
data/raw/methylation_summary.parquet
```

This table summarizes methylation by CpG site and tumor type, including tumor median, normal median, pan-cancer references, delta methylation, dispersion, and sample counts.

### 3. Gene and CpG annotation

```bash
python pipelines/06_build_gene_annotation.py
python pipelines/07_build_gene_map.py
python pipelines/08_build_cpg_gene_map.py
```

Expected outputs:

```text
data/raw/gene_annotation.parquet
data/raw/cpg_gene_map.parquet
```

These files define gene coordinates, TSS positions, and CpG-to-gene mappings.

### 4. Expression-methylation correlation

```bash
python pipelines/09_build_expr_parquet.py
python pipelines/10_cpg_expression_corr_tumor.py
python pipelines/11_merge_correlations.py
```

Expected output:

```text
data/raw/expression_correlation.parquet
```

This table contains CpG-level methylation-expression correlations per tumor type.

### 5. Additional biomarker features

```bash
python pipelines/12_build_cpg_feature.py
python pipelines/13_build_biomarker_candidates.py
python pipelines/14_generate_tumor_types.py
```

Expected outputs:

```text
data/raw/cpg_features.parquet
data/raw/tumor_types.parquet
```

`cpg_features.parquet` contains additional CpG-level information such as leukocyte methylation background.

### 6. PostgreSQL build

Set local PostgreSQL variables:

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

The build script creates the database if needed, drops existing project tables, recreates the schema, imports `.parquet` files, creates indexes, and runs `ANALYZE`.

## Notes

- The PostgreSQL build step is destructive for existing project tables because it drops and rebuilds them.
- Do not run `build_postgres.py` against a production database unless you intend to recreate the tables.
- Large input and output files should not be committed to GitHub.
