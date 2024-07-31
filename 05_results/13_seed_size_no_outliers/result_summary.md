# GWAS on seed size with and without a kinship matrix

**Date:** 17th July 2024
**Author:** Tom Ellis

## Background

Result 07_gemma_seed_size found a hint of a peak in F8 replicate 1, which we 
suspect is due to a small number of lines with very large seeds.

## What did you do?

- `02_library/run_GEMMA_with_without_K.sh` runs the GWAS
- `02_library/plot_gwas.py` plots it.
- `01_identify_outliers.sh` Plots distributions of seed size
- `02_runs_gemma.sh` runs the scripts on phenotype data.

## Main conclusion

Parents and rep1 have three outliers with BLUP values > 0.03
However, they aren't the same lines.
Removing the outliers removes the peak.

## Caveats

Not clear how these data were randomised.
F8s have no replication.
Parents and F8s not grown together.

## Follow-up

Improve resolution with haplotype blocks.