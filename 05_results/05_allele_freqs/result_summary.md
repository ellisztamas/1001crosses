# Changes in allele frequencies

**Date:** 8th July 2024
**Author:** Tom Ellis

## Background

Check that allele frequencies are largely the same between parents and F8s

## What did you do?

- `01_allele_freqs.sh` uses plink to get allele frequencies from VCF files.
- `02_plot_allele_freq_change.R` plots the above
- `03_af_from_gwas.R` pulls allele frequencies from a GWAS file

## Main conclusion

Raw allele frequencies are a mess.
They are well correlated when you use the GWAS file.
Seems like plink is doing something whacky.

## Caveats

Everything related to plink, to be honest

## Follow-up
