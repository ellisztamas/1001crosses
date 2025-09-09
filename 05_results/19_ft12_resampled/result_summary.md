# GWAS for flowering time on resampled data

**Date:** 17th June 2025
**Author:** Tom Ellis

## Background

Previous GWAS indicates that GWAS peaks are variable between parents and F9 cohorts.
That suggests GWAS is not very repeatable.
Since we have independent replicates of the crosses we can take independent
samples of these and repeat the GWAS for each.

## What did you do?

- Create 200 phenotype files by resampling from observed BLUPs
    - For lines that exist as a pair, take one of each
    - For singleton lines, take half
- Run GWAS on each.
- Merge files to keep only chr, ps, pvals and beta


## Main conclusion

Subsets tend to recapitulate the same peaks.

## Caveats


## Follow-up

