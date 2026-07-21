# Pipelines README — MethyMarker Biomarker Discovery App

This folder contains the reproducible data-processing workflow used by the methylation biomarker discovery app.

The workflow has two conceptual layers:

1. **CpG-level evidence**
   - tumor/normal methylation summaries;
   - pan-cancer methylation background;
   - PB methylation background;
   - CpG-to-gene mapping;
   - methylation-expression correlation.

2. **Region-level evidence**
   - tumor-independent physical CpG regions;
   - sequence-context features;
   - region-level prioritization in Streamlit using dynamic tumor-specific filters.

The region layer is intentionally separated into two steps: first define physical regions from the manifest, then add sequence features independently of tumor evidence. This avoids circularity and allows the same region universe to be reused across tumor types.

---

## Workflow overview

```text
Raw public methylation, phenotype, expression, PB and annotation data
        ↓
01 Clean CpG manifest
        ↓
02 Clean TCGA phenotype metadata
        ↓
03 Filter TCGA methylation matrix
        ↓
04 Compute tumor/normal and pan-cancer methylation summaries
        ↓
05 Merge cohort-level summaries
        ↓
06 Build gene annotation table
        ↓
07 Build gene symbol to Ensembl mapping
        ↓
08 Build CpG-to-gene mapping
        ↓
09 Preprocess RNA-seq expression matrix
        ↓
10 Compute CpG-gene methylation-expression correlations
        ↓
11 Merge cohort-level correlation files
        ↓
12 Build PB methylation background features
        ↓
13 Generate tumor-type lookup table
        ↓
14 Build tumor-independent physical CpG regions
        ↓
15 Add tumor-independent sequence features to regions
        ↓
build_postgres Load tables into PostgreSQL
        ↓
Streamlit app: Region Explorer and Gene Explorer
```

---

## Main outputs used by the app

| Output | Description |
|---|---|
| `data/raw/manifest_clean.parquet` | Clean CpG annotation table with genomic coordinates and gene annotation. |
| `data/raw/pheno_clean.parquet` | Clean TCGA phenotype metadata with tumor/normal sample classification. |
| `data/raw/methy.parquet` | Filtered TCGA methylation beta-value matrix. |
| `data/summary/{COHORT}_summary.parquet` | Cohort-specific CpG-level methylation summary statistics. |
| `data/raw/methylation_summary.parquet` | Unified tumor-summary table loaded as `tumor_summary`. |
| `data/raw/gene_annotation.parquet` | Gene annotation table with coordinates, strand, and TSS. |
| `data/raw/gene_map.parquet` | Gene symbol to Ensembl ID mapping table. |
| `data/raw/cpg_gene_map.parquet` | Expanded CpG-to-gene mapping table. |
| `data/raw/expr.parquet` | Filtered and log-transformed TCGA RNA-seq expression matrix. |
| `data/correlations_tumor/{COHORT}_Tumor_corr.parquet` | Cohort-specific CpG-gene Spearman correlations. |
| `data/raw/expression_correlation.parquet` | Unified methylation-expression correlation table. |
| `data/raw/cpg_features.parquet` | PB methylation reference features. |
| `data/raw/tumor_types.parquet` | Tumor-type lookup table. |
| `data/biomarkers/biomarker_region.parquet` | Tumor-independent physical CpG region table. |
| `data/biomarkers/biomarker_region_cpg.parquet` | Bridge table linking physical regions to CpGs. |
| `data/biomarkers/biomarker_region_sequence_score.parquet` | DNA sequence features and sequence-context score for each physical region. |
| `data/biomarkers/biomarker_region_sequence_score_config.json` | Parameters and diagnostics from sequence-score generation. |

---

## Script descriptions

### `generate_pb_methylation.R`

Generates PB methylation beta-value tables from raw GEO IDAT files using the `sesame` package.

**Input**

```text
data/geo/IDATS/
```

**Output**

A PB methylation matrix used by `12_build_cpg_feature.py`.

**Purpose**

This layer helps deprioritize CpGs or regions that are highly methylated in PB-derived DNA.

---

### `01_clean_manifiesto.py`

Downloads and preprocesses the Illumina HumanMethylation450K manifest.

**Processing**

- keeps relevant CpG annotation fields;
- retains valid autosomal CpG probes;
- removes probes without usable chromosome, position, or gene annotation;
- standardizes site IDs, gene symbols, CpG-island annotation, and genomic coordinates.

**Output**

```text
data/raw/manifest_clean.parquet
```

---

### `02_clean_phenotype.py`

Downloads and cleans TCGA phenotype metadata.

**Processing**

- keeps TCGA samples;
- retains primary tumor and solid tissue normal samples;
- generates `patient_id` from TCGA sample IDs;
- standardizes tumor type, tissue type, sample type, sex, age, and sample class.

**Output**

```text
data/raw/pheno_clean.parquet
```

---

### `03_methy_download_clean_v2.py`

Downloads and filters the TCGA methylation beta-value matrix.

**Inputs**

```text
data/raw/manifest_clean.parquet
data/raw/pheno_clean.parquet
```

**Processing**

- reads the methylation matrix in chunks;
- keeps CpGs present in the cleaned manifest;
- keeps samples present in the cleaned phenotype table;
- stores a compact Parquet matrix for downstream summary calculations.

**Output**

```text
data/raw/methy.parquet
```

---

### `04_methy_summary_v2.py`

Computes CpG-level methylation summary statistics by tumor type.

For each selected cohort, the script summarizes:

- tumor samples from the selected cohort;
- normal samples from the selected cohort;
- pan-cancer tumor samples excluding the selected cohort;
- pan-cancer normal samples excluding the selected cohort.

The main statistics are:

```text
tumor_median
tumor_std
tumor_n
normal_median
normal_std
normal_n
pan_tumor_median
pan_tumor_std
pan_tumor_n
pan_normal_median
pan_normal_std
pan_normal_n
delta_median
hi_index
```

**Output**

```text
data/summary/{COHORT}_summary.parquet
```

---

### `05_merge_methy_summary.py`

Merges cohort-specific methylation summaries into one unified table.

**Input**

```text
data/summary/*_summary.parquet
```

**Output**

```text
data/raw/methylation_summary.parquet
```

This table is loaded into PostgreSQL as `tumor_summary`.

---

### `06_build_gene_annotation.py`

Builds a gene annotation table with gene symbols, Ensembl IDs, genomic coordinates, strand, biotype, and TSS.

**Output**

```text
data/raw/gene_annotation.parquet
```

This table supports gene-centered visualization and coordinate-based sequence inspection.

---

### `07_build_gene_map.py`

Builds a gene symbol to Ensembl ID mapping table.

**Processing**

- extracts unique gene symbols from the manifest;
- queries gene ID resources;
- removes duplicated mappings;
- standardizes symbols and Ensembl IDs.

**Output**

```text
data/raw/gene_map.parquet
```

---

### `08_build_cpg_gene_map.py`

Builds the CpG-to-gene relationship table.

**Processing**

- expands CpGs annotated to multiple genes;
- generates one row per CpG-gene pair;
- adds Ensembl Gene IDs when available;
- standardizes gene symbols for exact joins in Streamlit.

**Output**

```text
data/raw/cpg_gene_map.parquet
```

Expected columns:

```text
site_id
gene_symbol
ensembl_id
```

---

### `09_build_expr_parquet.py`

Preprocesses TCGA RNA-seq expression data.

**Processing**

- keeps samples present in the cleaned phenotype table;
- keeps genes linked to CpGs in the methylation workflow;
- removes Ensembl version suffixes;
- applies `log2(FPKM-UQ + 1)` transformation;
- stores values as `float32` to reduce file size.

**Output**

```text
data/raw/expr.parquet
```

---

### `10_cpg_expression_corr_tumor.py`

Computes methylation-expression correlations in tumor samples.

**Processing**

- selects tumor samples for one cohort;
- matches methylation and expression values by shared TCGA sample IDs;
- computes Spearman correlation for each CpG-gene pair;
- stores correlation coefficient, p-value, and sample count.

**Output**

```text
data/correlations_tumor/{COHORT}_Tumor_corr.parquet
```

Expected columns:

```text
site_id
tumor_type
sample_class
ensembl_id
gene_symbol
spearman_r
pvalue
n_samples
```

If the cohort-specific output does not include `tumor_type` or `sample_class`, `11_merge_correlations.py` should add them during the merge step.

---

### `11_merge_correlations.py`

Merges cohort-specific methylation-expression correlation files.

**Input**

```text
data/correlations_tumor/*_corr.parquet
```

**Output**

```text
data/raw/expression_correlation.parquet
```

Important: this file should preserve the gene-specific fields:

```text
site_id
tumor_type
sample_class
ensembl_id
gene_symbol
spearman_r
pvalue
n_samples
```

Do not collapse only by `site_id` and `tumor_type`, because a CpG may map to more than one gene.

---

### `12_build_cpg_feature.py`

Computes PB methylation reference features.

**Processing**

- keeps CpGs present in the cleaned manifest;
- filters CpGs with excessive missingness;
- imputes remaining missing values using CpG-specific medians;
- computes PB methylation median and standard deviation.

**Output**

```text
data/raw/cpg_features.parquet
```

Expected columns:

```text
site_id
pb_median
pb_std
n_samples
```

---

### `13_generate_tumor_types.py`

Generates the tumor-type lookup table used by the app.

**Output**

```text
data/raw/tumor_types.parquet
```

Expected columns:

```text
tumor_type
full_name
tissue
```

---

### `14_build_biomarker_regions.py`

Builds tumor-independent physical CpG regions from the cleaned manifest.

This pipeline does **not** use tumor methylation, normal methylation, expression, PB features, or biological scores. It only uses CpG genomic coordinates and gene annotation. This separation is important because the same physical region universe can later be reused across tumor types without being biased by a specific cohort.

**Input**

```text
data/raw/manifest_clean.parquet
```

**Core logic**

1. Standardizes manifest columns to `site_id`, `gene_symbol`, `chr`, and `start_pos`.
2. Expands multi-gene annotations when needed.
3. Groups CpGs by gene and chromosome.
4. Creates clusters when consecutive CpGs are separated by no more than the selected maximum gap.
5. Requires a minimum number of CpGs per cluster.
6. Defines core region boundaries using the first and last CpG in the cluster.
7. Adds flanking browser coordinates.
8. Collapses duplicated physical intervals while preserving all associated genes.
9. Creates a bridge table linking each region to its CpGs.

**Outputs**

```text
data/biomarkers/biomarker_region.parquet
data/biomarkers/biomarker_region_cpg.parquet
data/biomarkers/biomarker_candidate_region_curve_config.json
```

**Useful columns in `biomarker_region.parquet`**

| Column | Meaning |
|---|---|
| `region_id` | Stable coordinate-derived physical region ID. |
| `gene_symbol` | Primary display gene. |
| `gene_symbols_all` | All gene symbols associated with the physical interval. |
| `chr` | Chromosome. |
| `core_start`, `core_end` | CpG-defined region boundaries. |
| `browser_start`, `browser_end` | Display interval with flanking bases. |
| `n_manifest_cpgs` | Number of CpG probes in the physical region. |
| `n_manifest_c` or `n_manifest_c_positions` | Number of unique CpG C-coordinate positions. |
| `cpg_density_per_100bp` | CpG probe density in the core region. |
| `region_rank_by_density` | Rank by CpG density. |

**Useful columns in `biomarker_region_cpg.parquet`**

| Column | Meaning |
|---|---|
| `region_id` | Physical region ID. |
| `site_id` | CpG probe ID. |
| `gene_symbol` | Gene symbol used for the region association. |
| `site_gene_symbols_all` | Full gene annotation for the CpG site. |
| `chr` | Chromosome. |
| `start_pos` | CpG genomic coordinate. |
| `cpg_order` | Order of the CpG within the region. |

**Default execution**

```bash
python pipelines/14_build_biomarker_regions.py
```

---

### `15_add_sequence_features_to_regions.py`

Adds DNA sequence features and a sequence-context score to physical regions generated by pipeline 14.

This pipeline is tumor-independent. It does not decide whether a region is a biomarker by itself. It provides a reusable sequence prior that can later be combined with tumor-specific methylation behavior, PB background, expression correlation, and region-level support.

**Inputs**

```text
data/biomarkers/biomarker_region.parquet
data/reference/hg38.fa
```

The FASTA must match the genome build used by the CpG manifest. For the current hg38 coordinates, use hg38.

**Core logic**

1. Loads physical regions from `biomarker_region.parquet`.
2. Uses either the core interval or browser interval depending on the configured sequence mode.
3. Reads the DNA sequence from a local FASTA.
4. Counts C, G, CG, and GCGC motifs.
5. Computes GC fraction, CG density, GCGC density, and manifest CpG density.
6. Normalizes density features using robust scaling.
7. Calculates `sequence_score`.
8. Writes a config JSON with execution parameters and diagnostics.

**Sequence score**

```text
sequence_score = 100 × (
    0.35 × normalized CG density
  + 0.30 × normalized GCGC density
  + 0.20 × GC fraction
  + 0.15 × normalized manifest CpG density
)
```

The score is not tumor-specific. It should be interpreted as a sequence-context prior.

**Outputs**

```text
data/biomarkers/biomarker_region_sequence_score.parquet
data/biomarkers/biomarker_region_sequence_score_config.json
```

**Useful columns**

| Column | Meaning |
|---|---|
| `region_id` | Physical region ID used for joins. |
| `gene_symbol` | Primary gene symbol, when preserved from the region table. |
| `chr` | Chromosome. |
| `sequence_region` | Whether core or browser coordinates were used. |
| `sequence_start`, `sequence_end` | Coordinates retrieved from FASTA. |
| `sequence_available` | Whether sequence retrieval succeeded. |
| `sequence_error` | Error message if sequence retrieval failed. |
| `sequence_length` | Length of retrieved sequence. |
| `n_c_sequence`, `n_g_sequence` | C and G counts. |
| `n_cg_sequence` | Overlapping CG count. |
| `n_gcgc` | Overlapping GCGC count. |
| `gc_fraction` | Fraction of bases that are C or G. |
| `cg_density_per_100bp` | CG motifs per 100 bp. |
| `gcgc_density_per_100bp` | GCGC motifs per 100 bp. |
| `manifest_cpg_density_per_100bp` | CpG probe density inherited from pipeline 14. |
| `sequence_score` | Final sequence-context score from 0 to 100. |

**Default execution**

```bash
python pipelines/15_add_sequence_features_to_regions.py --fasta data/reference/hg38.fa
```

---

### `build_postgres.py`

Creates the PostgreSQL schema and imports the processed Parquet tables.

**Expected tables**

```text
cpg_annotation
cpg_gene_map
expression_correlation
cpg_features
sample_metadata
gene_annotation
tumor_types
tumor_summary
biomarker_cpg_score
biomarker_region
biomarker_region_cpg
biomarker_region_sequence_score
```

**Local PostgreSQL**

```bash
export POSTGRES_USER=postgres
export POSTGRES_PASSWORD=your_password
export POSTGRES_HOST=localhost
export POSTGRES_PORT=5432
export POSTGRES_DB=db_methylation

python pipelines/build_postgres.py
```

**Remote PostgreSQL / Neon / other providers**

```bash
export DATABASE_URL="postgresql+psycopg2://user:password@host:5432/dbname"

python pipelines/build_postgres.py
```

---

## Recommended execution order

Run from the project root:

```bash
Rscript pipelines/generate_pb_methylation.R

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

---

## Region-level interpretation

The region-level workflow separates three ideas:

1. **Physical region definition**
   - where CpG-rich candidate intervals exist in the genome;
   - handled by pipeline 14;
   - tumor-independent.

2. **Sequence context**
   - whether the interval has favorable CG/GCGC/GC/CpG-density structure;
   - handled by pipeline 15;
   - tumor-independent.

3. **Tumor-specific evidence**
   - whether CpGs in the region are hypermethylated in tumor, low in normal tissues, low in PB, heterogeneous across tumor samples, and inversely correlated with expression;
   - handled dynamically in Streamlit.

This structure supports manuscript-level explanation: the genomic region universe is fixed first, sequence priors are added second, and biological/tumor evidence is applied afterward.

---

## Reproducibility notes

- Large raw and intermediate files should remain outside Git version control.
- Parquet outputs are used to reduce disk usage and improve downstream loading speed.
- Region IDs are derived from genomic coordinates and remain stable if input coordinates and parameters are unchanged.
- Multi-gene annotations are preserved rather than discarded.
- The sequence score depends on the selected FASTA and region mode; the FASTA path and parameters are stored in the config JSON.
- Pan-cancer reference summaries should exclude the selected cohort to avoid circular comparisons.
- Spearman correlation is used for methylation-expression association because the relationship is not assumed to be linear.
