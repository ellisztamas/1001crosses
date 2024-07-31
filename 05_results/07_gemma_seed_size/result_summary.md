# GWAS on seed size with and without a kinship matrix

**Date:** 27th June 2024
**Author:** Tom Ellis

## Background

A first attempt to run GWAS on seed size data in the parents and two replicate
cohorts of F9 plants, comparing results with and without adjusting for
population structure using a kinship matrix.

## What did you do?

- `02_library/run_GEMMA_with_without_K.sh` runs the GWAS
- `02_library/plot_gwas.py` plots it.
- `01_run_gemma.sh` runs the scripts on phenotype data.

## Main conclusion

Models without K matrix are inflated, especially the parents, as you would expect.

Using K matrix, there are no clear peaks in the parents or F8 rep2.
F8 rep1 shows a fairly broad peak on chr3, which I think is due to a haplotype
block in three outlier lines with especially large seeds, real or not.

## Caveats

Not clear how these data were randomised.
F8s have no replication.
Parents and F8s not grown together.

## Follow-up

Improve resolution with haplotype blocks.

Result 13_seed_size_no_outliers follows up on the peak in F8 rep1 excluding
the large lines.