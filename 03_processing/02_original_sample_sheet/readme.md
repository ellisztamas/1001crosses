Scripts to extract raw sequencing data for the F8s and call SNPs.

Script names should be fairly informative.
See the explanation within each script for more details.

Note that scripts 5, 6, 9, 10 and 11 currently don't work, and probably don't
need to, and that the sample names in the output VCF file are base file names.

To do:
    - Rename this something like "02_snp_calls"
    - Remove the SNPmatch stuff and move these to a folder for validating genotypes.
    - Either fix, or remove the base calibration stuff.
    - Add a script to rename the samples to give plate and positions.
        This will be enough to validate genotypes. Give the samples proper names
        only after validation.