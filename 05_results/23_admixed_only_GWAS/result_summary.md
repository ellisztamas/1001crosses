# GWAS for North-South intercross lines

**Date:** 13th October 2025
**Author:** Tom Ellis

## Background

Parents are mostly stratified by North vs South (PC1), then South 1 vs South 2
groups.
PCA for the F8s look like there is stratification both within along these axes.
Furthermore, GWAS for flowering time, seed size and rosette size don't return
anything convincing.
One explanation is that a lot of the crosses are *within* axdmiture groups, so
they don't really break up LD.
What happens if we only use lines that are crosses *between* admixture groups.

## What did you do?

- `01_subset_lines.qmd`: Uses the PC to identify lines descended from one
    Northern and one Southern parent.
- `02_run_gwas.sh`: Run GWAS for those lines only.
- `03_seed_size_conditional_GWAS.sh`: Repeat the GWAS, including the genotypes
    at the most strongly associated locus as a covariate.
- `04_seed_size_no_Groensvik.sh` GWAS excluding lines with the QTL haplotype.
- `05_extract_QTL_region.sh`: Script to extract SNP calls at the associated 
    locus as a matrix.
- `06_PCA_of_QTL_region.R`: PCA of the SNP calls in the associated region. 
- `07_condition_on_PC1.sh`: GWAS including PC1 of the QTL region as a covariate.    
- `08_kruskal_wallis_GWAS.sh`: Kruskall-Wallis GWAS using subset lines.
- `09_seed_size_GWAS.qmd` Summary of GWAS results.
- I attempted to look for structural variants, but this didn't work.

## Main conclusion

There is a 30kb haplotype on Chr5 associated with seed size.
The lines with that hapltoype seem to come from the Grönsvik peninsula.

If the association is real, this predicts that parental or other F9 lines with 
the haplotype should also have larger seeds.
However, they look like a random draw from the population.

Conditioning on the top SNP returns a slightly larger haplotype.
Removing the Grönsvik lines returns nothing convincing.
Kruskall Wallis GWAS returns nothing.

## Caveats

I'm not 100% sure GEMMA is including covariates in the way I think it is,
although it reports that it is.

## Follow-up

Structural variants.