#' Prepare phenotype data for GWAS.
#'
#' Import data on number of days to flowering at 10°C for the F8s, and results
#' from SNPmatch confirming that these are the genotypes we think they are.
#' For correct genotypes, create separate files for replicates 1 and 2.
#'
#' Input: Flowering-time data, SNPmatch results
#' Output: PLINK formatted phenotype and covariate files for GEMMA
#'
#' Tom Ellis, 5th December 2023.
library(tidyverse)

# Directory for saving the output files (which will be the input for GEMMA)
outdir <- "05_results/04_gwas_ft10/gemma_input_files"
dir.create(outdir, showWarnings = FALSE)

# Import data
# Import SNPmatch results
sm_results <- read_csv(
  "../crosses_continued/004.F8/002.mapping/001.data/plate_info.csv",
  show_col_types = FALSE
) %>%
  rename(
    sample = Sample_name,
    CorrectParents = genoSelection
    ) %>%
  filter(
    CorrectParents == "yes"
  ) %>%
  select(sample, CorrectParents)
# Phenotype data
# Merge with SNPmatch results, filter for correct matches and split the sample name
ft10 <- read_csv(
  "01_data/05_phenotype_expt/flowering_time_F8s.csv",
  col_types = 'cfffcci'
  ) %>%
  mutate(
    sample = gsub(' ', '_', sample) # Change spaces to underscores in the variable `name`
  ) %>%
  right_join( sm_results, by = 'sample') %>%
  separate(sample, into=c("cross",'generation','rep'), remove=FALSE) %>%
  filter( !is.na(FT10) )

# Write phenotype files for reps 1 and 2 to disk
# This is in PLINK format, which requires columns for familyID, within-familyID and phenotype.
# We don't have families, so the first column is a column of zeroes
ft10 %>%
  filter(rep == "rep1") %>%
  mutate(dummy_column = 0) %>%
  select(dummy_column, sample, FT10) %>%
  write_tsv(paste0(outdir, "/ft10_rep1.tsv"), col_names = FALSE)
ft10 %>%
  filter(rep == "rep2") %>%
  mutate(dummy_column = 0) %>%
  select(dummy_column, sample, FT10) %>%
  write_tsv(paste0(outdir, "/ft10_rep2.tsv"), col_names = FALSE)

#' Write covariate files for reps 1 and 2 to disk
#' The first column is a column of ones to tell Gemma to fit an intercept
#' The second column experimental block (1 to 5)
ft10 %>%
  filter(rep == "rep1") %>%
  mutate(dummy_column = 1) %>%
  select(dummy_column, plate) %>%
  write_tsv(paste0(outdir, "/ft10_rep1_covariates.tsv"), col_names = FALSE)
ft10 %>%
  filter(rep == "rep2") %>%
  mutate(dummy_column = 1) %>%
  select(dummy_column, plate) %>%
  write_tsv(paste0(outdir, "/ft10_rep2_covariates.tsv"), col_names = FALSE)
