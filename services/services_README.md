# Services

This folder contains external service utilities used by the app.

## `ensembl.py`

Utilities for retrieving genomic sequences from Ensembl or related genome APIs.

Used by the Gene Explorer sequence browser to display genomic context around selected CpG sites.

## Design note

External API calls should be cached when possible to reduce latency and avoid repeated requests for the same genomic region.
