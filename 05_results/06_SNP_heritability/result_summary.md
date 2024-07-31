# SNP heritability

**Date:** 4th July 2024
**Author:** Tom Ellis

## Background

Calculate SNP heritability for rosette size and seed size based on BLUPs.
See 03_processing/06_process_phenotypes for how these are estimated

## What did you do?

- `02_library/run_SNP_heritability.sh ` is a general script to create a relatedness
    matrix and calculate the variance explained using GEMMA.
- `01_calculate_h2.sh` runs the script.

## Main conclusion

Heritabilities are 0.2 or lower.
For both traits, F8 rep1 was lowest, F8 rep2 highest, with parents in the middle.

## Caveats

Seed size data are external

## Follow-up

Repeat with seed size and flowering time data from my phenotyping experiment when
these become available.