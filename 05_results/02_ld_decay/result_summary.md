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

## Main conclusion

Short-range LD looks identical between the F8s and parents - this is good

## Caveats

## Follow-up

Work out a sensible way to plot long-range LD.