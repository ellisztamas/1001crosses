# Validation of the F8 genotypes

This directory contains attempts to validate the F8 genotypes by comparing them to known genotypes of the parents in the published SNP matrix.

1. Using SNPmatch
2. Visually, using ibdpainting.

The most informative results are summarised in `06_ibd_painting_results.md` and the HTML it generates.

Things to do next:
- Regmap
- Directory for the new sequencing data
    Gives a VCF file with 91 samples of F8s
    Another with two parents
- Directory for breakpoints
    - Swap names for the six samples
    - Run rTiger on those samples
        - Work out whether I need to merge 1137 and 1074 with the other parents
    - Slice up the VCF file one sample at a time