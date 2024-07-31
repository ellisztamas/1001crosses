# GWAS on rosette size conditioned on top SNP

**Date:** 23rd July 2024
**Author:** Tom Ellis

## Background

result 08 suggested there's a peak on chr3 for rosette size.
Repeat the GWAS including the genotype at that SNP as a cofactor to see if any
other peaks pop out

## What did you do?

- `02_library/run_GEMMA_with_without_K.sh` runs the GWAS
- `02_library/plot_gwas.py` plots it.
- `01_extract_top_SNP.py` extracts the genotype at the top SNP and saves as a text file.
- `02_run_gemma.sh` runs the scripts on phenotype data.

## Main conclusion

The peak is still there. This is very strange.

## Caveats

May be that nearby SNPs are apparently not correlated because of SNP-call errors

## Follow-up

Repeat this after hapotype calls.