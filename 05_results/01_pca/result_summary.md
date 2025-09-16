# Population structure in parents and crosses

**Date:** 18th December 2023
**Author:** Tom Ellis

## Background

How has the structure changed between the parents and F8s?

## What did you do?

- `01_run_pca.sh` prepares and runs a PCA on VCF files for the parents and F8s.
- `02_plot_structure.R` plots the first 2 eigenvectors
- `03_pca_genetic_variation.qmd` and the resulting HTML file investigate the PCs
    in more detail, and provide correlates with other variables.

## Main conclusion

- Parents are stratified along two axes.
    - PC1 is largely correlated with latitude.
    - This is not perfect, check metadata
- F8s are are dispersed along those axes, and between.
- Flowering time somewhat correlated with a bunch of PCs.
- Seed size shows substantial correlation with both PC1 and PC2
  - Even in the F8s!
  - Suggests seed size is highly polygenic?

## Caveats

## Follow-up

- Resequence 1435, 6199 and 5835.
