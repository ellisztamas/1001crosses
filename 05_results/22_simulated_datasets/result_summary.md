# Simulated GWAS datasets uses GCTA

**Date:** 30th July 2025
**Author:** Tom Ellis

## Background

What causes GWAS statistics to be inflated?

- Long-range LD
- Short-range LD
- Number of loci
- Heritability?
- Selection? (not run)
- Epistasis (not run)

I want to simulate a bunch of GWAS phenotypes, run GWAS, and calculate the 
lambda statistic on each.

## What did you do?

- `01_permute_VCF.py` is a CLI script to randomise genotypes at each row in a
    VCF file to completely remove linkage disequilibrium, even at short scales. 
- `02_permute_parents.sh` applies the Python script to the parental VCF file.
- `03_create_plink.sh` converts VCF files for observed (non-permuted) VCF files
    for the parents and F8s to Plink format.
- `04_simulate_phenotypes.sh` uses GCTA to simulate phenotypes for each dataset.
- `05_run_gemma.sh` runs GWAS on the simulated datasets.
- `06_calculate_lambda.R` is a CLI script to load each GEMMA results file
    and calculate the lambda statistic and the number of 500kb windows with at
    least one Bonferroni significant association.
- `07_calculate_lambda.sh` runs the Rscript on the cluster.


## Main conclusion

Heritability is the main driver of inflation, as measured by lambda.
Polygenicity has an effect but it is remarkably weak.

LD also has a huge effect.
Parents show massive inflation, but also high variance in lambda, often being
below 1.

## Caveats

This was done before I updated the imputation. The scripts won't work out of
the box because they use separate VCF files for F8 cohorts 1 and 2.

It looks like QQ plots are often pathological, so I don't know if one can really
use lambda as a good summary statistic.

## Follow-up

Maybe I should just look at lambda -log10 p < 2.
