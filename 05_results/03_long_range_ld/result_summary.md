# Plot the decay in LD

**Date:** 28th April 2025
**Author:** Tom Ellis

## Background

Plot long-range linkage disequilibrium (LD) between loci in the genome.

The goal is to create a heatmap of the LD matrix like that in Long et al. (2013)
for only pairs of SNPs that are above a certain threshold (Long used 0.6).
The heatmap is masked to show only the values above the threshold.
The upper triangle should show values for the parents, and the lower triangle
values for the progeny.

The challenge is that for 10e6 SNPs, the LD matrix is 10e12 values.

## What did you do?

- `01_sequential_pruning.sh` uses PLINK to filter SNPs to look at. Based on a 
    suggestion from Yoav Voichek, I used the `--indep-pairwise` to remove SNPs
    in *local* LD with one another (within 1kb). This is a greedy algorithm
    that goes from SNP to SNP, identifies pairs in LD, and keeps only the marker
    with the higher MAF. This saves a list of SNPs for later.
- `02_LD_between_windows.py` is a library script that takes those SNP positions,
    a thousand at a time, and looks up their genotypes in the imputed HDF5
    genotype files. It calculates LD between pairs of these chunks, and stores
    positions and r2 between pairs of SNPs with LD > a threshold.
- `03_find_SNPs_in_LD.sh` executes the Python script on the data. The array is
    over chunks of 1000 SNPs in the pruned list.
- `04_plot_heatmap.R` Plots pairwise LD for parents and progeny.
    Parents are shown in the upper left, progeny in the lower right.

## Main conclusion

As expected, there is tons of LD>0.5 between chromosomes in the parents
There is essentially none in the progeny.
What there is is r2 < 0.3.

## Caveats

- I've used the list of SNPs retained from the parents only.
    One could think about a better way to do this.
    The logic was that LD decays within 250 kb in A. thaliana, so 1kb seems 
    kinda fine.
- This is with imputed data - is that appropriate?

## Follow-up