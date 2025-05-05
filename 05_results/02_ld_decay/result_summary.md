# Plot the decay in LD over short distances

**Date:** 5th May 2024
**Author:** Tom Ellis

## Background

Linkage disequlibrium after one round of random mating ought to be half that of
parents. Compare the decay of LD over short distances.

## What did you do?

- `01_get_LD.sh` calculates LD in plink.
- `02_plot_LD.R` plots the results

## Main conclusion

LD starts lower and decays faster in the progeny than in the parents.
At distances of 100kb LD levels out at about
    r2=0.07 in the parents
    r2=0.05 in the progeny
So there is less of a difference than one would expect.
But LD in the progeny seems to dip and come back up at around 30kb

## Caveats

LD in the progeny seems to dip, and also starts lower than I would expect.

## Follow-up