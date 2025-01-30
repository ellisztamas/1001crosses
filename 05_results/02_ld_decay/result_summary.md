# Plot the decay in LD

**Date:** 27th June 2024
**Author:** Tom Ellis

## Background

Linkage disequlibrium after one round of random mating ought to be half that of
parents. Check that this is true by comparing F8s and parents.

## What did you do?

- `01_get_LD.sh` calculates LD in plink2.
- `02_plot_LD.R` plots the results
- `03_genome_wide_LD.sh` an attempt to assess long range LD. This didn't work,
    because the output is massive.
- Two python scripts calculate LD directly using HDF5 files for 10000 evenly
    space SNPs

## Main conclusion

Short-range LD looks identical between the F8s and parents - this is good
Most Long-range LD is on the order of 0.1 to 0.3 (weak) for these loci
There are fewer SNPs in long-range LD for F9 cohort 1 than for the parents

## Caveats

## Follow-up

Need to revisit this when I've got SNP matrices in order.
For long-range LD it would be better to first go through each pair of loci and
    identify those with LD above some threshold (Long eta al used 0.6), then
    *only plot these*.
I will have to roll my own way to do this, because PLINK gets OOM killed.