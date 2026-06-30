# MethyMarker — Cancer Methylation Biomarker Discovery App

Interactive and reproducible platform for prioritizing DNA methylation biomarker candidates across TCGA tumor cohorts. The project integrates tumor/normal methylation summaries, pan-cancer background, leukocyte methylation reference profiles, CpG-to-gene annotation, methylation-expression correlation, physical CpG regions, and sequence-context scores into a PostgreSQL-backed Streamlit application.

The current workflow is **region-first**. CpG sites remain the primary evidence unit, but the user-facing biomarker candidate is a **physical CpG region** supported by one or more qualifying CpGs.

---

## Project overview

DNA methylation is frequently altered in cancer and can provide candidate biomarkers when tumor specificity, normal-tissue background, leukocyte background, gene context, expression association, and local sequence structure are considered together.

This repository contains:

- reproducible preprocessing pipelines;
- derived Parquet tables;
- a PostgreSQL schema/loading script;
- Streamlit pages for region-level prioritization and gene-level inspection.

The application supports:

- tumor-vs-normal methylation comparison;
- pan-cancer background comparison;
- leukocyte methylation filtering;
- CpG-to-gene mapping;
- methylation-expression correlation;
- physical CpG-region aggregation;
- sequence-context scoring using CG/GCGC/GC features;
- interactive candidate prioritization for paper-oriented biomarker discovery.

---

## Repository structure

```text
biomarker_discovery_app/
│
├── Main.py
│
├── pages/
│   ├── 1_Region_Explorer.py
│   ├── 2_Gene_Explorer.py
│   └── pages_README.md
│
├── pipelines/
│   ├── generate_leukocytes_methylation.R
│   ├── 01_clean_manifiesto.py
│   ├── 02_clean_phenotype.py
│   ├── 03_methy_download_clean_v2.py
│   ├── 04_methy_summary_v2.py
│   ├── 05_merge_methy_summary.py
│   ├── 06_build_gene_annotation.py
│   ├── 07_build_gene_map.py
│   ├── 08_build_cpg_gene_map.py
│   ├── 09_build_expr_parquet.py
│   ├── 10_cpg_expression_corr_tumor.py
│   ├── 11_merge_correlations.py
│   ├── 12_build_cpg_feature.py
│   ├── 13_generate_tumor_types.py
│   ├── 14_build_biomarker_regions.py
│   ├── 15_add_sequence_features_to_regions.py
│   ├── build_postgres.py
│   └── pipelines_README.md
│
├── db/
│   ├── db.py
│   ├── queries.py
│   └── db_README.md
│
├── data/
│   ├── raw/
│   ├── summary/
│   ├── correlations_tumor/
│   ├── biomarkers/
│   ├── reference/
│   └── geo/
│       └── IDATS/
│
├── services/
│   └── ensembl.py
│
├── requirements.txt
├── environment.yml
├── env.example
├── schema.sql
├── schema_neon.sql
└── README.md
```

Large raw and intermediate data files are not expected to be versioned in Git. They should be downloaded or generated locally and stored under `data/`.

---

## Data sources

| Data layer | Source | Role in the app |
|---|---|---|
| TCGA phenotype metadata | UCSC Xena / GDC Pan-Cancer | Defines tumor type, tissue class, sample type, and patient/sample IDs. |
| TCGA methylation | UCSC Xena / GDC Pan-Cancer HumanMethylation450K | Provides CpG beta-values used for tumor/normal and pan-cancer summaries. |
| TCGA expression | UCSC Xena / GDC Pan-Cancer RNA-seq | Used to estimate methylation-expression association for CpG-gene pairs. |
| CpG annotation | HumanMethylation450K manifest mapped to hg38 | Provides CpG IDs, genomic coordinates, gene annotations, and CpG island context. |
| Leukocyte methylation | GEO IDAT files processed with `sesame` | Provides a blood-cell background filter for biomarker specificity. |
| Gene annotation | GENCODE / Ensembl-derived annotation | Provides gene coordinates, strand, TSS, and gene symbols. |
| DNA sequence | Local hg38 FASTA | Used to compute sequence-context features for physical CpG regions. |

---

## Conceptual workflow

```text
Raw public datasets and annotation
        ↓
Clean CpG manifest and TCGA phenotype metadata
        ↓
Filter TCGA methylation beta-value matrix
        ↓
Compute tumor/normal and pan-cancer CpG summaries
        ↓
Build gene annotation and CpG-to-gene maps
        ↓
Preprocess expression and compute CpG-gene correlations
        ↓
Generate leukocyte methylation background features
        ↓
Build tumor-independent physical CpG regions
        ↓
Add tumor-independent sequence-context features
        ↓
Load normalized tables into PostgreSQL
        ↓
Prioritize candidate regions in Streamlit
        ↓
Inspect selected genes, regions, CpGs, and sequence context
```

The separation between **physical region definition**, **sequence-context scoring**, and **tumor-specific evidence** is intentional. It avoids using tumor behavior to define the genomic universe and makes the model easier to explain in a manuscript.

---

## Main derived outputs

| Output | Description |
|---|---|
| `data/raw/manifest_clean.parquet` | Clean CpG manifest with site IDs, genomic coordinates, and gene annotation. |
| `data/raw/pheno_clean.parquet` | Clean TCGA phenotype metadata with tumor/normal classification. |
| `data/raw/methy.parquet` | Filtered methylation beta-value matrix. |
| `data/summary/{COHORT}_summary.parquet` | Cohort-level CpG methylation summary statistics. |
| `data/raw/methylation_summary.parquet` | Unified tumor-summary table loaded as `tumor_summary`. |
| `data/raw/gene_annotation.parquet` | Gene-level annotation table. |
| `data/raw/gene_map.parquet` | Gene symbol to Ensembl ID mapping. |
| `data/raw/cpg_gene_map.parquet` | Expanded CpG-to-gene bridge table. |
| `data/raw/expr.parquet` | Filtered and log-transformed expression matrix. |
| `data/raw/expression_correlation.parquet` | CpG-gene methylation-expression correlation table. |
| `data/raw/cpg_features.parquet` | Leukocyte methylation background features. |
| `data/raw/tumor_types.parquet` | Tumor-type lookup table. |
| `data/biomarkers/biomarker_region.parquet` | Tumor-independent physical CpG region table. |
| `data/biomarkers/biomarker_region_cpg.parquet` | Bridge table linking regions to CpGs. |
| `data/biomarkers/biomarker_region_sequence_score.parquet` | Sequence features and sequence score for each region. |
| `data/biomarkers/biomarker_region_sequence_score_config.json` | Parameters used for sequence scoring. |

---

## PostgreSQL model

The Streamlit app expects the processed Parquet outputs to be loaded into PostgreSQL.

### Base tables

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
biomarker_region
biomarker_region_cpg
biomarker_region_sequence_score
biomarker_cpg_score
```

The database loader is:

```bash
python pipelines/build_postgres.py
```

The loader supports two connection modes:

```bash
# Local PostgreSQL
export POSTGRES_USER=postgres
export POSTGRES_PASSWORD=your_password
export POSTGRES_HOST=localhost
export POSTGRES_PORT=5432
export POSTGRES_DB=db_methylation

python pipelines/build_postgres.py
```

or:

```bash
# Direct database URL
export DATABASE_URL="postgresql+psycopg2://user:password@host:5432/dbname"

python pipelines/build_postgres.py
```

---

## Running the preprocessing workflow

Run scripts from the project root.

```bash
Rscript pipelines/generate_leukocytes_methylation.R

python pipelines/01_clean_manifiesto.py
python pipelines/02_clean_phenotype.py
python pipelines/03_methy_download_clean_v2.py
python pipelines/04_methy_summary_v2.py
python pipelines/05_merge_methy_summary.py
python pipelines/06_build_gene_annotation.py
python pipelines/07_build_gene_map.py
python pipelines/08_build_cpg_gene_map.py
python pipelines/09_build_expr_parquet.py
python pipelines/10_cpg_expression_corr_tumor.py
python pipelines/11_merge_correlations.py
python pipelines/12_build_cpg_feature.py
python pipelines/13_generate_tumor_types.py
python pipelines/14_build_biomarker_regions.py
python pipelines/15_add_sequence_features_to_regions.py --fasta data/reference/hg38.fa
python pipelines/build_postgres.py
```

The FASTA used in pipeline 15 must match the genome build of the CpG coordinates. For the current hg38 manifest, use an hg38 FASTA.

---

## Running the Streamlit app

Configure the database URL either in `.streamlit/secrets.toml`:

```toml
[database]
url = "postgresql+psycopg2://user:password@host:5432/dbname"
```

or through the environment:

```bash
export DATABASE_URL="postgresql+psycopg2://user:password@host:5432/dbname"
```

Then run:

```bash
streamlit run Main.py
```

---

## Streamlit pages

### `Main.py`

Landing page for the application. It points users to the active modules:

```text
Region Explorer
Gene Explorer
```

### `pages/1_Region_Explorer.py`

Global region-level prioritization page. It applies tumor-specific CpG filters and aggregates qualifying CpGs into physical regions. It produces bubble plots where each point is a region.

Main score:

```text
sequence_site_score = sequence_score × n_qualifying_sites
```

When expression is enabled:

```text
expression_signal = max(0, -mean_spearman_r)
expression_score_component = 100 × expression_signal × n_qualifying_sites
final_region_score = sequence_site_score + expression_score_component
```

When expression is disabled:

```text
final_region_score = sequence_site_score
```

### `pages/2_Gene_Explorer.py`

Gene-level validation page. It shows the methylation profile across the selected gene, overlays filtered candidate regions, lists region/CpG evidence, and opens a sequence browser for selected regions or sites.

---

## Reproducibility notes

- Large raw and derived data files should remain outside Git version control.
- Parquet outputs are used for compact storage and fast loading.
- Methylation beta-values are summarized using medians to reduce sensitivity to outliers.
- Pan-cancer reference groups should exclude the cohort being analyzed.
- Spearman correlation is used for methylation-expression association because the relationship is not assumed to be linear.
- Region IDs are derived from genomic coordinates and are stable when input coordinates and parameters are unchanged.
- Sequence score depends on the FASTA and region mode used by pipeline 15; these parameters are stored in the config JSON.

---

## Suggested methods statement

TCGA HumanMethylation450K beta-values, phenotype metadata, and RNA-seq expression data were processed together with CpG manifest annotation, leukocyte methylation reference profiles, and gene annotation. CpG-level tumor/normal methylation summaries were computed by tumor type and compared against pan-cancer background summaries. CpG-gene methylation-expression associations were estimated using Spearman correlation in tumor samples. Tumor-independent physical CpG regions were generated from manifest coordinates and annotated with sequence-context features from an hg38 reference FASTA. Processed outputs were loaded into PostgreSQL and explored through an interactive Streamlit application for region-level biomarker prioritization and gene-level validation.

---

## Repository status

This repository is intended for research and reproducibility. The current version supports exploratory and paper-oriented prioritization of methylation biomarker candidate regions across selected TCGA tumor cohorts.
