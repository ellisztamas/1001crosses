#' Prepare flowering-time data the parents for GWAS.
#'
#' Import data on number of days to flowering at 10°C for the parents and format
#' for running GEMMA.
#'
#' To do: consider subsetting the parents to match those in each offspring file.
#' This may not be necessary, because if we split F8s by replicate sample sizes
#' are very similar.
#'
#' Input: Flowering-time data for the F8s and parents
#' Output: PLINK formatted phenotype and covariate files for GEMMA
#'
#' Tom Ellis, 12th December 2023.

library(tidyverse)

# Directory for saving the output files (which will be the input for GEMMA)
outdir <- "05_results/04_gwas_ft10/gemma_input_files"
dir.create(outdir, showWarnings = FALSE)

# Flowering time data for the parents
parents <- read_csv(
  "01_data/05_phenotype_expt/flowering_time_10C_parents.csv",
  col_types="cd"
  )

# List of parents in the F8s
# I had thought to subset the data to match, but they are similar enough not to
# worry about this for now.
read_tsv(
  "05_results/04_gwas_ft10/gemma_input_files/ft10_rep2.tsv",
  col_names = c("null_column","sample", "FT10C"),
  col_types = 'icd'
  ) %>%
  select(sample, FT10C) %>%
  separate(sample, sep="_", into=c("cross", "generation", "rep")) %>%
  separate(cross, into=c("parent1", 'parent2'), sep="x") %>%
  pivot_longer(parent1:parent2) %>%
  pull(value) %>%
  unique()

parents %>%
  filter(!is.na(FT_10C)) %>%
  mutate(dummy_column = 0) %>%
  select(dummy_column, accession, FT_10C) %>%
  write_tsv(paste0(outdir, "/ft10_parents.tsv"), col_names = FALSE)
