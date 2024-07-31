# Plan for checking genotypes in a transparent way

If we are to use these data we need to validate identities in a transparent way.

Why can there be mismatches?

- Samples are mixed up in the 1163 SNP matrix
- Material is mixed up:
    - There were mix-ups in the lines that were crossed
    - Seed stocks have become confused in the subsequent 8 generations
    - Samples were mixed/contaminated at sequencing
- There are dubious SNPs in the 1163 matrix, meaning genotype calls are whacky
- Sequencing data for the F8s are poor
- Processing or genotype calls for the F8s are dubious

How to address that:

1. Be conservative about which SNPs to use for genotype calls.
2. Make a table of the lines for whom the most likely parents are not the expected parents
    - Count up how often each appears
3. Compare that list to the list of 'mixed-up' lines from Rahul's paper.
4. Check whether any line is consistently coming up as something else
    - If yes, this was probably a mix-up at crossing
    - Do a visual check of similarity across these genomes
    - Cross check with Rahul's mix-ups
5. Check whether the expected parent comes up within top candidates, even if not first.
    - Visual check for where the differences are (there may to too many for this to be realistic)
    - Set some criteria for how different they can be
    - Quantify relateness between the expected and most likely parent
        - Prediction is they should be closely related
6. Verify all this by comparing F8s to F9s
    - (If we ran the collection through PhenoPlant we'd need to sequence F9s anyway)