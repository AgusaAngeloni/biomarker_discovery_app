# Pipeline README — Biomarker Discovery App

This document summarizes the data-processing workflow used by the methylation biomarker discovery app. The repository is organized around a reproducible pipeline that starts from public methylation, phenotype, expression, leukocyte and annotation resources, generates normalized Parquet tables, loads them into PostgreSQL, and exposes the results through Streamlit visualizations.

The current pipeline has two conceptual layers:

1. **CpG-level evidence**: tumor/normal methylation summaries, pan-cancer background, leukocyte methylation background and methylation-expression correlation.
2. **Region-level evidence**: physical CpG-rich regions built from the manifest, sequence-context features, and downstream prioritization of regions that contain CpGs passing biologically meaningful filters.

The region layer is intentionally separated into two steps: first define physical regions independently of tumor evidence, then add sequence features independently of tumor evidence. This avoids circularity and allows the same region universe to be reused across tumor types.

---

## Workflow overview

```text
Raw public data and annotation
        ↓
Clean CpG manifest and phenotype metadata
        ↓
Filter TCGA methylation matrix
        ↓
Compute tumor/normal and pan-cancer methylation summaries
        ↓
Build gene and CpG-to-gene mapping tables
        ↓
Preprocess RNA-seq expression matrix
        ↓
Compute CpG-expression correlations in tumor samples
        ↓
Compute leukocyte methylation background features
        ↓
Build tumor-independent physical CpG regions
        ↓
Add tumor-independent sequence features to regions
        ↓
Load processed tables into PostgreSQL
        ↓
Explore CpGs, genes and candidate regions in Streamlit
```

---

## Main outputs used by the app

| Output | Description |
|---|---|
| `data/raw/manifest_clean.parquet` | Clean CpG annotation table with genomic coordinates and gene annotation. |
| `data/raw/pheno_clean.parquet` | Clean TCGA phenotype metadata with tumor/normal sample classification. |
| `data/raw/methy.parquet` | Filtered TCGA methylation beta-value matrix. |
| `data/summary/{COHORT}_summary.parquet` | Cohort-specific CpG-level methylation summary statistics. |
| `data/raw/methylation_summary.parquet` | Unified tumor-summary table across cohorts. |
| `data/raw/gene_map.parquet` | Gene symbol to Ensembl ID mapping. |
| `data/raw/cpg_gene_map.parquet` | Expanded CpG-to-gene mapping table. |
| `data/raw/expr.parquet` | Filtered and log-transformed expression matrix. |
| `data/raw/expression_correlation.parquet` | CpG-gene Spearman correlations in tumor samples. |
| `data/raw/cpg_features.parquet` | Leukocyte methylation reference features. |
| `data/biomarkers/biomarker_region.parquet` | Tumor-independent physical CpG region table. |
| `data/biomarkers/biomarker_region_cpg.parquet` | Bridge table linking physical regions to CpGs. |
| `data/biomarkers/biomarker_region_sequence_score.parquet` | DNA sequence features and sequence score for each physical region. |
| `data/raw/tumor_types.parquet` | Tumor-type lookup table used by the app. |

---

## Pipeline scripts

### `generate_leukocytes_methylation.R`

Generates leukocyte methylation beta-value tables from raw GEO IDAT files using the `sesame` package. This creates the leukocyte reference layer used to avoid prioritizing CpGs that are highly methylated in blood-derived cells.

**Input**

- Raw IDAT files under `data/geo/IDATS/`.

**Output**

- Leukocyte methylation matrix used by `12_build_cpg_feature.py`.

---

### `01_clean_manifiesto.py`

Downloads and preprocesses the Illumina methylation manifest.

**Processing**

- Keeps relevant CpG annotation fields.
- Retains valid CpG probes.
- Removes probes without usable chromosome, position or gene annotation.
- Standardizes site IDs, gene symbols and genomic coordinates.

**Output**

- `data/raw/manifest_clean.parquet`

---

### `02_clean_phenotype.py`

Downloads and cleans TCGA phenotype metadata.

**Processing**

- Keeps TCGA samples.
- Retains primary tumor and solid tissue normal samples.
- Generates `patient_id` from sample IDs.
- Standardizes tumor type, tissue type and sample class.

**Output**

- `data/raw/pheno_clean.parquet`

---

### `03_methy_download_clean_v2.py`

Downloads and filters the TCGA methylation beta-value matrix.

**Inputs**

- TCGA methylation matrix.
- `manifest_clean.parquet`.
- `pheno_clean.parquet`.

**Processing**

- Reads the methylation matrix in chunks.
- Keeps CpGs present in the cleaned manifest.
- Keeps samples present in the cleaned phenotype table.
- Stores a compact Parquet matrix for downstream summary calculations.

**Output**

- `data/raw/methy.parquet`

---

### `04_methy_summary_v2.py`

Computes CpG-level methylation summary statistics by tumor type.

**Processing**

For each tumor cohort, the script summarizes beta-values for:

- tumor samples from the selected cohort;
- normal samples from the selected cohort;
- pan-cancer tumor samples excluding the selected cohort;
- pan-cancer normal samples excluding the selected cohort.

The main statistics include median beta-value, standard deviation and valid sample count. Tumor/normal delta and dispersion metrics are then used by the app and biomarker-prioritization steps.

**Output**

- `data/summary/{COHORT}_summary.parquet`

---

### `05_merge_methy_summary.py`

Merges cohort-specific methylation summaries into a unified table.

**Output**

- `data/raw/methylation_summary.parquet`

---

### `07_build_gene_map.py`

Builds a gene identifier conversion table.

**Processing**

- Extracts unique gene symbols from the manifest.
- Retrieves Ensembl Gene IDs.
- Removes duplicated mappings.

**Output**

- `data/raw/gene_map.parquet`

---

### `08_build_cpg_gene_map.py`

Builds the CpG-to-gene relationship table.

**Processing**

- Expands CpGs annotated to multiple genes.
- Generates one row per CpG-gene pair.
- Adds Ensembl Gene IDs when available.

**Output**

- `data/raw/cpg_gene_map.parquet`

---

### `09_build_expr_parquet.py`

Preprocesses TCGA RNA-seq expression data.

**Processing**

- Keeps samples present in the cleaned phenotype table.
- Keeps genes linked to CpGs in the methylation workflow.
- Removes Ensembl version suffixes.
- Applies `log2(FPKM-UQ + 1)` transformation.
- Stores values as `float32` to reduce file size.

**Output**

- `data/raw/expr.parquet`

---

### `10_cpg_expression_corr_tumor.py`

Computes methylation-expression correlations in tumor samples.

**Processing**

- Selects tumor samples for one cohort.
- Matches methylation and expression by shared TCGA sample IDs.
- Computes Spearman correlation for each CpG-gene pair.
- Stores correlation coefficient, p-value and sample count.

**Output**

- `data/correlations_tumor/{COHORT}_Tumor_corr.parquet`

---

### `11_merge_correlations.py`

Merges cohort-specific methylation-expression correlation files.

**Output**

- `data/raw/expression_correlation.parquet`

---

### `12_build_cpg_feature.py`

Computes leukocyte methylation reference features.

**Processing**

- Keeps CpGs present in the cleaned manifest.
- Filters CpGs with excessive missingness.
- Imputes remaining missing values using CpG-specific medians.
- Computes leukocyte methylation median and standard deviation.

**Output**

- `data/raw/cpg_features.parquet`

---

### `13_generate_tumor_types.py`

Generates the tumor type lookup table used by the app.

**Output**

- `data/raw/tumor_types.parquet`

---

### `14_build_biomarker_regions_2.py`

Builds tumor-independent physical CpG regions from the cleaned manifest. This is the first region-level step of the biomarker workflow.

This pipeline does **not** use tumor methylation, normal methylation, expression, leukocyte features or biological scores. It only uses CpG genomic coordinates and gene annotation. This separation is important because the same physical region universe can later be reused across tumor types without being biased by a specific cohort.

**Input**

- `data/raw/manifest_clean.parquet`

**Core logic**

1. Standardizes manifest columns to `site_id`, `gene_symbol`, `chr` and `start_pos`.
2. Expands multi-gene annotations by default, so a CpG annotated to `GENE1;GENE2` can contribute to both gene-specific cluster searches.
3. Groups CpGs by gene and chromosome.
4. Creates candidate clusters when consecutive CpGs are separated by no more than `--max-gap-bp` and the cluster contains at least `--min-cpgs` CpGs.
5. Builds core region coordinates from the first and last CpG in the cluster.
6. Adds flanking browser coordinates using `--flank-bp`.
7. Collapses duplicated physical intervals by genomic coordinates, preserving all associated genes in `gene_symbols_all`.
8. Creates a normalized bridge table linking regions to CpG sites.

**Outputs**

- `data/biomarkers/biomarker_region.parquet`
- `data/biomarkers/biomarker_region_cpg.parquet`

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
| `n_manifest_c_positions` | Number of unique CpG genomic positions. |
| `cpg_density_per_100bp` | CpG probe density in the core region. |
| `region_rank_by_density` | Rank by CpG density. |

**Default execution**

```bash
python pipelines/14_build_biomarker_regions_2.py
```

**Common options**

```bash
python pipelines/14_build_biomarker_regions_2.py   --max-gap-bp 350   --min-cpgs 2   --flank-bp 100   --gene-mode expand
```

---

### `15_add_sequence_features_to_regions_v2.py`

Adds DNA sequence features and a sequence-context score to the physical regions generated by pipeline 14.

This pipeline is also tumor-independent. It does not decide whether a region is a biomarker by itself. Instead, it provides a reusable sequence prior that can later be combined with tumor-specific methylation behavior, leukocyte background, expression correlation and biological filters.

**Input**

- `data/biomarkers/biomarker_region.parquet`
- Local uncompressed reference FASTA passed with `--fasta`

**Important genome-build note**

The FASTA must match the coordinate system used by the manifest. In this case use hg38.

**Core logic**

1. Loads physical regions from `biomarker_region.parquet`.
2. Uses either the core interval or the browser interval depending on `--sequence-region`.
3. Reads the sequence from a local FASTA.
4. Counts C, G, CG and GCGC motifs.
5. Computes GC fraction, CG density, GCGC density and manifest CpG density.
6. Normalizes density features using robust 95th-percentile scaling.
7. Calculates `sequence_score`.
8. Writes a config JSON with the execution parameters and sequence diagnostics.

**Sequence score**

```text
sequence_score = 100 × (
    0.35 × normalized CG density
  + 0.30 × normalized GCGC density
  + 0.20 × GC fraction
  + 0.15 × normalized manifest CpG density
)
```

The score is not tumor-specific. It should be interpreted as a sequence-context prior for biomarker prioritization.

**Output**

- `data/biomarkers/biomarker_region_sequence_score.parquet`
- `data/biomarkers/biomarker_region_sequence_score_config.json`

**Useful columns**

| Column | Meaning |
|---|---|
| `region_id` | Physical region ID used for joins. |
| `sequence_region` | Whether core or browser coordinates were used. |
| `sequence_start`, `sequence_end` | Coordinates retrieved from FASTA. |
| `sequence_available` | Whether sequence retrieval succeeded. |
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
python pipelines/15_add_sequence_features_to_regions_v2.py   --fasta data/reference/hg19.fa
```

**Core-only sequence score**

```bash
python pipelines/15_add_sequence_features_to_regions_v2.py   --fasta data/reference/hg19.fa   --sequence-region core
```

**Browser-window sequence score**

```bash
python pipelines/15_add_sequence_features_to_regions_v2.py   --fasta data/reference/hg19.fa   --sequence-region browser
```

---

## Recommended execution order

Run from the project root:

```bash
Rscript pipelines/generate_leukocytes_methylation.R

python pipelines/01_clean_manifiesto.py
python pipelines/02_clean_phenotype.py
python pipelines/03_methy_download_clean_v2.py
python pipelines/04_methy_summary_v2.py
python pipelines/05_merge_methy_summary.py
python pipelines/07_build_gene_map.py
python pipelines/08_build_cpg_gene_map.py
python pipelines/09_build_expr_parquet.py
python pipelines/10_cpg_expression_corr_tumor.py
python pipelines/11_merge_correlations.py
python pipelines/12_build_cpg_feature.py
python pipelines/13_generate_tumor_types.py
python pipelines/14_build_biomarker_regions_2.py
python pipelines/15_add_sequence_features_to_regions_v2.py --fasta data/reference/hg19.fa
```

Then build or refresh the PostgreSQL database using the database-loading scripts configured for the project.

---

## Region-level interpretation

The region-level workflow separates three ideas that should not be mixed too early:

1. **Physical region definition**: where CpG-rich candidate intervals exist in the genome. This is handled by pipeline 14.
2. **Sequence context**: whether the interval is enriched in potentially informative sequence features such as CG/GCGC content. This is handled by pipeline 15.
3. **Tumor-specific evidence**: whether CpGs in the region are hypermethylated in tumor, low in normal tissues, low in leukocytes, heterogeneous across tumor samples, and inversely correlated with expression. This is handled by downstream scoring and the Streamlit pages.

This structure keeps the workflow reproducible and makes it easier to explain in a manuscript: the genomic region universe is fixed first, sequence priors are added second, and biological/tumor evidence is applied afterward.

---

## Reproducibility notes

- All large raw and intermediate files should remain outside Git version control.
- Parquet outputs are used to reduce disk usage and improve downstream loading speed.
- Region IDs are derived from genomic coordinates, making them stable across runs if the input coordinates and parameters are unchanged.
- Multi-gene annotations are preserved rather than discarded.
- The sequence score depends on the selected FASTA and genome build; the FASTA path and mode are stored in the config JSON.
- Pan-cancer reference summaries should exclude the selected cohort to avoid circular comparisons.
- Spearman correlation is used for methylation-expression association because the relationship is not assumed to be linear.
