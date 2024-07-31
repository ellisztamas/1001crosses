# Simulate GWAS datasets with a single causative SNP

**Date:** Early July 2024
**Author:** Tom Ellis

## Background

In the absence of decent phenotype data I want to generate simulations.
The idea is to generate datasets where we know the correct answer examine what
the K-matrix correction does, and whether it does what we think.

- Use the observed VCF files for parents and F8s
- Choose a SNP to be a real QTL
- Simulate a binary phenotype based on genotypes at that QTL with probabilities
    0.5, 0.6, 0.7 and 0.8
- Run GEMMA
- Quantify the p-value, estimated effect size and genomic inflation factor

## What did you do?

I previously picked SNPs to simulate in result 09_gwas_single_loci; see scripts 
01 to 03 in that result.

- `01_simulate_phenotypes.py` creates phenotype files in plink format
- `02_submit_simulate_phenotypes.sh` runs the above as a SLURM job.
- `02_library/plot_gwas.py` plots it.
- `03_run_gemma.sh` runs GEMMA on simulated datasets with and without the
    correction for the relatedness matrix.
- `04_summarise_simulations.sh` calls `02_library/summarise_sims.py` on each 
    GWAS result and creates a text file with a single line.
- `05_concatenate_summaries.sh` concatenates all the single-line text files into
    one usable CSV file.

## Main conclusion

Uncorrected GWAS are inflated, even when there is no background effect.
The K matrix largely fixes this.

Associations are stronger for stronger QTL (logical)

## Caveats

## Follow-up

Compare with simulations that include background effects (results 11 and 12).