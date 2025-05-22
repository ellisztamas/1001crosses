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
- F8s are more homogeneous, except for 6 samples that have whacky parents.
- Flowering time correlated with a bunch of PCs.

## Caveats

- PC1 should probably be more strongly correlated tbh
- Offspring of parents from the RegMap panel make a cluster of their own

## Follow-up

- Check the GPS metadata.
- Resequence 1435, 6199 and 5835.