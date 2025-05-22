# GWAS for flowering time at 12°C

**Date:** 22nd May 2025
**Author:** Tom Ellis

## Background

GWAS with GEMMA on BLUPs for log flowering time.

## What did you do?

- `01_run_gemma.sh` runs GEMMA with and without a K matrix
- `02_plot_gwas.Rmd` plots all Manhattan and QQ plots
- `03_plot_parents_vs_cohorts.R` makes some basic plots comparing pvalues and
    effect sizes in the parents and F8s
- `04_gwas_statistics.qmd` is a more detailed examination of GWAS statistics for
    models with no relatdness matrix.

## Main conclusion

Without using a relatedness matrix:

- Inflation in the F9s is about half that of the parents (lambda ~ 2 vs 4).
- That's still really high! Probably polygenicity.
- Effect sizes strongly correlated, pvals less so.
- Effects in parents show some outliers
- Results don't look very repeatable

## Caveats


## Follow-up

Plot candidate genes.

Why do pvals in the parents show outliers?
Check allele frequencies or genetic variances
    Are they the outliers in the bottom left of the AF change plot?
Plot them on the manhattan plot - same place?

Take subsets of F9s and see how repeatable peaks are

Kruskall-Wallis GWAS