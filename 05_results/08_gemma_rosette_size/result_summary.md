# GWAS on rosette size with and without a kinship matrix

**Date:** 27th June 2024
**Author:** Tom Ellis

## Background

A first attempt to run GWAS on rosette size data in the parents and two replicate
cohorts of F9 plants, comparing results with and without adjusting for
population structure using a kinship matrix.

## What did you do?

- `02_library/run_GEMMA_with_without_K.sh` runs the GWAS
- `02_library/plot_gwas.py` plots it.
- `01_run_gemma.sh` runs the scripts on phenotype data.

## Main conclusion

Without K matrix, the parents are hugely inflated.
F8s are inflated, but less so, with clearer peaks.

Including K fixes the inflation, possibly overcorrecting.
The parents and F8 rep2 show no peaks.
F8 rep1 shows a wide peak on chr3, which I suspect is something low frequency.

## Caveats

## Follow-up

Condition on the top SNP in chromosome 3; see result 14_rosette_size_condition_on_top_SNP

Refine haplotypes.