# Streamlit Pages README

This folder contains the interactive Streamlit pages for the MethyMarker biomarker discovery app.

The current application logic is **region-first**:

```text
CpG-level tumor evidence
        ↓
Physical CpG region model
        ↓
Sequence-context score
        ↓
Interactive region prioritization
        ↓
Gene-level validation and sequence inspection
```

CpG sites remain the biological evidence unit, but the final user-facing candidate object is a **physical CpG region**.

---

## Active pages

```text
pages/
├── 1_Region_Explorer.py
├── 2_Gene_Explorer.py
└── pages_README.md
```

---

## Required PostgreSQL tables

The pages are designed to query PostgreSQL dynamically instead of loading whole Parquet files.

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

### Region / biomarker tables

```text
biomarker_region
biomarker_region_cpg
biomarker_region_sequence_score
biomarker_cpg_score
```

`tumor_summary` should expose the standardized names:

```text
hi_index
pan_tumor_median
pan_normal_median
```

Compatibility aliases such as `dispersion_index`, `HI_index`, `pantumor_median`, `pannormal_median`, `panTumor_median`, and `panNormal_median` may be handled by the page code, but the preferred schema uses the standardized names above.

---

## `1_Region_Explorer.py`

### Purpose

`1_Region_Explorer.py` is the global candidate-prioritization page. It starts from tumor-specific CpG filters and aggregates qualifying CpGs into physical candidate regions.

Each plotted point represents:

```text
one physical CpG region
```

not one CpG site.

This page is intended to answer:

- Which regions have multiple tumor-specific CpGs?
- Which regions are low in normal tissue and PB?
- Which regions have high tumor heterogeneity or dispersion?
- Which regions have favorable sequence context?
- Which regions show inverse methylation-expression association?

---

## Main sidebar filters

Typical filters include:

```text
Tumor type
Minimum delta_median
Maximum normal_median
Maximum pan_normal_median
Maximum pb_median
Minimum hi_index
Optional expression filter using mean Spearman r
```

The filters are applied at the CpG-evidence level. The passing CpGs are then aggregated into regions.

---

## Region aggregation logic

The page joins:

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
mean_pb_median
mean_spearman_r
sequence_score
cpg_sites
```

The page preserves the original `biomarker_region.region_id` as the stable region identifier. Coordinate-derived labels are only auxiliary labels.

---

## Region score logic

The base score favors regions with both strong sequence context and multiple qualifying CpGs:

```text
sequence_site_score = sequence_score × n_qualifying_sites
```

Expression is optional.

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

Interpretation:

- `sequence_score` is tumor-independent and comes from pipeline 15.
- `n_qualifying_sites` is tumor/filter-dependent and comes from the selected Streamlit thresholds.
- `mean_spearman_r` is tumor-specific and should ideally be negative for genes where methylation is associated with lower expression.

---

## Main plots

The page shows three region-level bubble plots:

1. **DMR plot**
   - X axis: mean HI index in the region
   - Y axis: mean delta methylation
   - Bubble size: number of qualifying CpGs

2. **Expression plot**
   - X axis: mean HI index in the region
   - Y axis: mean delta methylation
   - Bubble size: expression signal

3. **Final score plot**
   - X axis: mean HI index in the region
   - Y axis: mean delta methylation
   - Bubble size: final region score
   - Color: final region score

All plots are configured for SVG export through the Plotly toolbar.

---

## Region candidate table

The region table should include the fields needed for manuscript-oriented inspection:

```text
gene_region_id
gene_main
chr
browser_start
browser_end
final_region_score
sequence_score
n_qualifying_sites
n_manifest_cpgs
fraction_qualifying_sites
mean_delta
mean_hi
mean_normal_median
mean_pan_normal_median
mean_pb_median
mean_spearman_r
expression_signal
sequence_site_score
cpg_sites
region_ids
```

`cpg_sites` allows users to inspect which CpGs support the candidate region.

---

## `2_Gene_Explorer.py`

### Purpose

`2_Gene_Explorer.py` is the gene-level validation and sequence-inspection page.

It is intended to answer:

- What is the full methylation profile of a selected gene?
- Which physical regions overlap that gene?
- Which candidate regions pass the selected tumor-specific filters?
- Which CpGs support those regions?
- What sequence context surrounds the selected region or site?

---

## Gene Explorer behavior

For a selected tumor type and gene, the page:

1. loads all CpGs mapped to the selected gene;
2. plots the complete methylation profile by genomic coordinate;
3. overlays candidate physical regions;
4. applies region-candidate filters dynamically;
5. lists candidate regions and candidate CpGs;
6. opens a sequence browser for a selected region or selected CpG site.

The full gene methylation profile remains visible. Candidate regions are highlighted on top of that profile rather than replacing it.

---

## Gene-level methylation profile

The profile typically includes:

```text
Tumor median methylation
Normal median methylation
Pan-cancer tumor median methylation
Pan-cancer normal median methylation
PB median methylation
Delta methylation
HI index
Spearman expression correlation
```

Candidate regions are shown as a visual track or highlighted interval on the same coordinate system.

---

## Candidate region filters in Gene Explorer

The region filters mirror the logic of the Region Explorer:

```text
Candidate minimum delta_median
Candidate maximum normal_median
Candidate maximum pan_normal_median
Candidate maximum pb_median
Candidate minimum hi_index
Optional maximum Spearman r
```

The filtered regions are calculated live from PostgreSQL. No precomputed candidate-region curve table is required.

---

## Sequence browser

The sequence browser is region-aware and can display:

```text
browser window: browser_start to browser_end
core region: core_start to core_end
custom flank: user-defined flank around the selected region or site
```

The browser highlights:

- manifest CpGs;
- the selected CpG or candidate region;
- GCGC motifs;
- nucleotide identity;
- genomic coordinates.

The sequence browser uses `services/ensembl.py` for sequence retrieval in the Streamlit page. Pipeline 15 uses a local FASTA to calculate the region sequence score.

---

## Design notes

- PostgreSQL is the source of truth for the pages.
- Pages should apply SQL filters before loading results into pandas.
- CpG-level evidence should be aggregated into region-level summaries for prioritization.
- Physical `region_id` values from `biomarker_region` should be preserved.
- Coordinate-derived IDs should only be auxiliary display labels.
- Region visualizations should avoid duplicated physical regions caused by multi-gene annotations.
- Sequence scores come from `biomarker_region_sequence_score`, generated by pipeline 15.
- The final user-facing candidate object is a region, not an isolated CpG.
