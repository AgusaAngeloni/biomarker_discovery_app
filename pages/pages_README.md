# Streamlit pages

This folder contains the interactive pages of the MethyMarker app.

## `1_CpG_Explorer.py`

CpG-level exploration page.

Typical features:

- Tumor type selection
- Delta methylation filters
- Normal methylation filters
- Leukocyte methylation filters
- Expression correlation filters
- Candidate CpG ranking
- Interactive plots and downloadable result tables

## `2_Gene_Explorer.py`

Gene-level exploration page.

Typical features:

- Gene search
- CpG sites mapped to selected gene
- Methylation curves across CpGs
- Selection of CpG sites from plots or tables
- Genomic sequence browser around selected CpG sites

## Design note

Pages should query PostgreSQL dynamically using filters selected by the user. Avoid loading full database tables into pandas before filtering.
