# Jointly estimate all SNP effects at once

**Date:** 4th November 2025
**Author:** Tom Ellis

## Background

I want to convey a sense for how frequent false positive associations should be
for different levels of confounding with LD.
I previously did this for many levels of polygenicity (result 22), but this can
be challenging to get your head around.
I also want to check how MAF affects things, and that is a bit much with 
different levels of polygenicity.
Here I will do some more detailed simulations, but for the case where there is
only a single causal locus.

## What did you do?

- `01_sample_genotypes.sh` stratifies VCF files into minor-allele-frequency
    bins, and takes a random sample of markers from each.
- `02_simulate_phenotypes.R` is a CLI script to subset those markers and
    simulate phenotypes via a linear model. The same model is used for parents
    and F8s.
- `03_simulate_phenotypes.sh` runs the R script on the cluster.
- `04_run_gemma.sh` runs GWAS on the simulated datasets.
- `05_calculate_inflation.R` is a CLI script to load each GEMMA results file
    and calculate the lambda statistic and the number of 500kb windows with at
    least one Bonferroni significant association.
`06_calculate_inflation.sh` runs the Rscript on the cluster.


## Main conclusion

## Caveats

## Follow-up

