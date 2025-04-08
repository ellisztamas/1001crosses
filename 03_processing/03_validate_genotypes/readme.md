# Validation of the F8 genotypes

This directory contains attempts to validate the F8 genotypes by comparing them
to known genotypes of the parents in the published SNP matrix.

1. Using SNPmatch
2. Visually, using ibdpainting.

The most informative results are summarised in `06_ibd_painting_results.md` and
the HTML it generates.

`07_new_sample_named.R` and `08_correct_split_VCFs.sh` rename the samples in the
VCF files from plate positions (e.g. 1_G4) to the validated genotype names.

Scripts 09 and 10 run ibdpainting again for the ambiguous samples, but this time
against the RegMap panel.

Script 11 subsets the parental genotypes to include only valid parents of at
one line, and only those SNPs present in the F8s.