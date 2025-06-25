# Kruskall-Wallis GWAS

**Date:** 20th June 2025
**Author:** Tom Ellis

## Background

I previously ran parametric GWAS with and without a relatedness matrix.
Parametric analyses do Things to the results, so here I try a non-parametric
analysis.

## What did you do?

- `02_library/kruskall_wallis_GWAS.py` is a python script to loop over every
    SNP and perform a Kruskall-Wallis test at each.
- `./01_run_KW_GWAS.sh` runs the script on phenotypes.
- `./02_parametric_vs_nonparametric.qmd` and the HTML it renders compare results
    with parametric results from GEMMA.

## Main conclusion

Results from Kruskall-Wallis tests look very similar to parametric tests when
you don't use a K matrix.

## Caveats

I didn't compare with parametric results with a K matrix, but the comparison
is essentially the same as between GEMMA results with and without K.

## Follow-up

