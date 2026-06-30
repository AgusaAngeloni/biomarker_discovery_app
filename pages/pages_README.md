# Streamlit pages

This folder contains the interactive pages of the MethyMarker app.

The current app logic is **region-first**. CpG sites are still the biological evidence unit, but the main candidate objects shown to the user are now **physical CpG regions** generated from the manifest and linked to tumor-specific CpG evidence through PostgreSQL.

## Current page workflow

```text
CpG-level methylation evidence
        ↓
Physical region model
        ↓
Sequence context score
        ↓
Interactive region prioritization
        ↓
Gene-level inspection and sequence browser
```

The pages should not load whole database tables into memory. They should query PostgreSQL dynamically using the filters selected by the user.

## Required PostgreSQL tables

The pages assume that the normalized database model has already been loaded.

### Base tables

```text
cpg_annotation
cpg_gene_map
tumor_summary
cpg_features
expression_correlation
gene_annotation
sample_metadata
tumor_types
```

### Biomarker / region tables

```text
biomarker_cpg_score
biomarker_region
biomarker_region_cpg
biomarker_region_sequence_score
```

`tumor_summary` should use the standardized column name:

```text
hi_index
```

not `HI_index`.

---

## `1_Region_Explorer.py`

Region-level candidate exploration page.

### Purpose

This page starts from CpG-level tumor filters and aggregates the matching CpGs into physical candidate regions. Each plotted point represents **one physical region**, not one CpG site.

The page is intended for global region prioritization across genes and tumor types.

### Main filters

Typical filters include:

- Tumor type
- Minimum `delta_median`
- Maximum `normal_median`
- Maximum pan-cancer normal methylation
- Maximum leukocyte methylation
- Minimum `hi_index`
- Optional expression filter based on negative Spearman correlation

### Region aggregation logic

The page joins CpG evidence with the normalized region model:

```text
tumor_summary
    + cpg_annotation
    + cpg_gene_map
    + cpg_features
    + expression_correlation
    + biomarker_region_cpg
    + biomarker_region
    + biomarker_region_sequence_score
```

For each physical region, it summarizes:

```text
n_qualifying_sites
n_manifest_cpgs
fraction_qualifying_sites
mean_delta
mean_hi
mean_normal_median
mean_pan_normal_median
mean_leukocyte_median
mean_spearman_r
sequence_score
```

### Score logic

The current region score is designed to favor regions with both:

1. multiple qualifying CpGs, and
2. strong sequence context from pipeline 18.

```text
sequence_site_score = sequence_score × n_qualifying_sites
```

Expression is optional. When enabled:

```text
expression_signal = max(0, -mean_spearman_r)
expression_score_component = 100 × expression_signal × n_qualifying_sites
final_region_score = sequence_site_score + expression_score_component
```

When expression is disabled:

```text
final_region_score = sequence_site_score
```

### Main plots

The page shows three region-level bubble plots:

1. **Bubble size = number of qualifying CpGs**  
   Useful to identify regions supported by multiple tumor-specific CpGs.

2. **Bubble size = expression signal**  
   Useful to identify regions where methylation is inversely associated with expression.

3. **Bubble size/color = final region score**  
   Useful for region prioritization.

All plots are exported as SVG through the Plotly toolbar.

### Output table

The region candidate table includes the main biological, sequence, and score fields. It also includes `cpg_sites`, allowing the user to inspect which CpGs support each region.

---

## `2_Gene_Explorer.py`

Gene-level validation page with filtered candidate regions and sequence browser.

### Purpose

This page focuses on one selected gene and tumor type. It keeps the gene methylation profile as the main visualization, but overlays only the **filtered physical candidate regions** that pass the selected thresholds.

It does **not** require a precomputed `biomarker_candidate_region_curve` table or parquet file.

### Main behavior

For a selected gene, the page:

1. loads all CpGs mapped to the gene,
2. plots the methylation profile across genomic position,
3. applies candidate-region filters dynamically,
4. marks filtered candidate regions on the methylation profile,
5. lists candidate regions and candidate CpGs,
6. opens a sequence browser for the selected candidate region.

### Candidate region filters

The sidebar includes region-candidate filters such as:

```text
Candidate minimum Delta
Candidate max Median NT
Candidate max PanCancer NT
Candidate max Leukocytes
Candidate min HI
Optional candidate expression filter
```

The filtered regions are calculated live from PostgreSQL using the normalized region tables.

### Methylation profile

The main profile includes:

```text
Median Type T
Median Type NT
Median PanCan T
Median PanCan NT
Median Leukocytes
Delta
Optional HI
Expression profile plot
```

Candidate regions are shown as a blue track or dotted boundaries on top of the full-gene methylation profile. The curves are not restricted to the selected region; the full gene profile remains visible.

### Region-level table

The region table summarizes candidate regions for the selected gene and tumor type.

Typical fields include:

```text
region_id
region_gene_symbol
chr
core_start
core_end
browser_start
browser_end
n_qualifying_cpgs
n_manifest_cpgs
fraction_qualifying_cpgs
sequence_score
mean_delta
max_delta
mean_hi
max_hi
mean_tumor_median
mean_normal_median
mean_leukocyte_median
mean_spearman_r
qualifying_cpg_sites
qualifying_cpg_positions
```

Selecting one row opens that region in the sequence browser.

### Sequence browser

The browser is region-based rather than CpG-centered.

Available windows:

```text
browser       → browser_start/browser_end
core          → core_start/core_end
custom flank  → user-defined flank around the core region
```

The browser highlights:

- manifest CpGs in the selected region,
- the core region,
- GCGC/CGCG motifs,
- genomic coordinates,
- downloadable SVG sequence view.

The browser is intended for sequence-level inspection of the candidate region after methylation-based filtering.

---

## Optional / legacy pages

### `3_Region_Universe_Browser.py`

This page is useful for debugging or inspecting the full region universe produced by pipeline 16. It reads region-universe style data and opens a region sequence browser without requiring tumor-specific CpG filters.

Use it when checking:

- whether regions were generated correctly,
- whether `browser_start/browser_end` are appropriate,
- whether manifest CpGs are correctly located in the sequence,
- whether GCGC/CGCG motifs are detected as expected.

For the main analysis workflow, prefer:

```text
1_Region_Explorer.py
2_Gene_Explorer.py
```

---

## Recommended analysis order

```text
1. Run preprocessing and load PostgreSQL.
2. Run pipeline 14 to generate CpG biological scores.
3. Run pipeline 16 to generate physical CpG regions.
4. Run pipeline 18 to generate sequence scores per region.
5. Rebuild/load PostgreSQL biomarker tables.
6. Use 1_Region_Explorer.py to prioritize regions globally.
7. Use 2_Gene_Explorer.py to inspect selected genes and candidate regions.
```

---

## Design notes

- PostgreSQL is the source of truth for the pages.
- Pages should use SQL filters before loading data into pandas.
- CpG-level evidence should be aggregated into region-level summaries for prioritization.
- Physical region IDs from `biomarker_region.region_id` should be preserved.
- Coordinate-derived IDs may be used only as auxiliary labels.
- Region visualizations should avoid duplicating regions caused by multiple gene annotations.
- Sequence scores come from `biomarker_region_sequence_score`, generated by pipeline 18.
- The final user-facing candidate object is a **region**, not an isolated CpG.
