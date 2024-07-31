# Population structure in parents and crosses

**Date:** 18th December 2023
**Author:** Tom Ellis

## Background

How has the structure changed between the parents and F8s?

## What did you do?

- `01_run_pca.sh` prepares and runs a PCA on VCF files for the parents and both 
    replicates of the F8s using PLINK.
- `02_plot_structure.R` plots the first 2 eigenvectors

## Main conclusion

- Parents are stratified along two axes.
- Offspring show a similar shape, but with scatter perpendicular to the axes and
    a certain amount of intermediates.
- More scatter in PC3

## Caveats

## Follow-up
