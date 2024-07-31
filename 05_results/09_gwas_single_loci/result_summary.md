# Simulated GWAS datasets using SNPs as traits

**Date:** 8th July 2024
**Author:** Tom Ellis

## Background

I want to investigate how true/false positive associations are affected by 
population structure, and the correction for it.

In this simulation I will use genotypes at individual SNPs as phenotypes.

I want to:
    - select SNPs that are associated/not associated with population structure
    - do this for SNPs of different frequencies
    - Extract 
    - Run GWAS with and without the kinship matrix
        - Repeat in the F8s.
    - Quantify:
        - p-value for the causal SNP
        - Number of false positives
            - Related: genomic inflation factor
        - Estimated effect size of the causal SNP

## What did you do?

- `01_allele_correlations.py` creates a TSV file giving information to choose 
    which SNPs to use. See the doc string of that script for what that means.
- `02_submit_allele_correlations.sh` runs the above as a SLURM job.
- `03_choose_SNPs.R` creates something like a BED file giving positions of SNPs to test.
- `04_simulate_phenotypes.py` creates GEMMA-formatted phenotypes from SNP loci.
- `05_submit_simulate_phenotypes.sh`: Script to run 04_simulate_phenotypes.py as a SLURM job
- `06_run_GEMMA.sh`: Run GEMMA on simulated datasets with and without the
    correction for the relatedness matrix.
- `07_summarise_sims.py`: A poorly documented script to extract the information
    from GWAS results.
    
## Main conclusion

This did not work satisfactorily.
GEMMA tended to fail, which I think is because the phenotypes were too tightly
correlated with relatedness matrix (I saw PVE~1).

## Caveats

## Follow-up

Here, effect sizes are determined by allele frequencies. Try to vary effect sizes.

This should also fix the crashes