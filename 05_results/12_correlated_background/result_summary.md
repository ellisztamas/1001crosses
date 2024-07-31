# Simulate GWAS datasets with historical selection

**Date:** 17th July 2024
**Author:** Tom Ellis

## Background

Follow up on the simulations in results 10 and 11, but including a 
genetic background effect on top of the main SNP which is correlated with
population structure.
I choose 10000 random SNPs and give each a small effect on the (logit) phenotype.
The direction of effects at background effects is correlated with the first PC
of population sturcture in the parents, corresponding to the case where there
selection has acted.

See also 05_results/10_sim_case_controls/result_summary.md

## What did you do?

Scripts follow the same pattern as 05_results/10_sim_case_controls; see the
result summary in that directory.

The main difference is inclusion of background effects in 01_simulate_phenotypes.py.

## Main conclusion

There's more inflation than for the single-QTL model

Associations at main QTL are weaker, and effect sizes slightly underestimates
than in the single-QTL model, suggesting the GWAS is overcorrecting.

## Caveats

This is a very crude way of cobbling a background effect together.
Be stricter about the variance explained by the main QTL and background effects.

## Follow-up
