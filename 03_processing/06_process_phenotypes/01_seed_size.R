#' Prepare raw seed size data for further analysis.
#'
#' Estimate genetic values for each line from the phenotype data by fitting
#' BLUPs.
#'
#' Inputs:
#'    Raw data file for seed size (with one row per seed) including replicate lines 1 and 2.
#'      These are F9 seeds, and would be siblings of the plants grown in the chamber experiment
#' Outputs:
#'    Tab-separated text files formatted for GEMMA.
#'    The first column is a vector of 1s to tell GEMMA to fit an intercept.
#'    Subsequent columns show genotype and BLUP estimates of
#'    genetic values of seed size. There are separate file for replicate lines
#'    1 and 2, plus the parents.
#'
#' Tom Ellis, 25th June 2024

library(tidyverse)
library(lme4)

source("03_processing/04_pieters_VCF/04_align_original_corrected_sample_names.R")

# Import size data for each seed separately
seed_size_F9 <- read_delim("01_data/05_phenotypes/raw_data_seed_size_F9.csv", delim = ";")
# Tidy the plantID names up
seed_size_F9 <- seed_size_F9 %>%
  rename(genotype = `Pick description`) %>%
  filter(
    ! grepl( 'Aa1', genotype ), # Remove entries for A. arenosa
    genotype != "6025x6046  MIX of replica 1 and 2" # Boxeed mixed one pair of replicates. Remove.
  ) %>%
  mutate(
    # Swap the one or two spaces in genotypes for underscores.
    genotype = str_replace(genotype, "  ", "_"),
    genotype = str_replace(genotype, " ", "_")
  )

# Correct sample names and remove dubious genotypes.
seed_size_F9 <- seed_size_F9 %>%
  inner_join(sample_name_corrections, by=c("genotype" = "Sample_name_original")) %>%
  rename(
    genotype_original = genotype,
    genotype = Sample_name_corrected
  )

# Vector of unique accession IDs among all the parents for which we have data.
# This will be used to subset the 1001 genomes data
parents_to_include <- seed_size_F9 %>%
  separate(genotype, into=c("genotype", 'replicate')) %>%
  separate(genotype, into=c("parent1", 'parent2'), sep="x") %>%
  pivot_longer(parent1:parent2) %>%
  pull(value) %>%
  unique()

# Imprort seed size data for the 1001 genomes accessions, and subset to those
# that are parents of the crosses
seed_size_parents <- read_delim(
  # "01_data/05_phenotypes/1001_2014_2017_2022_boxeed.csv",
  '/groups/nordborg/projects/seed_size/01_data/boxeed/1001_2014_2017_2022_boxeed.csv',
  delim = ";"
) %>%
  rename(genotype = `Pick description`) %>%
  filter(
    genotype %in% parents_to_include
  )


# Fit BLUPs
seedsize_F9_blups <- ranef(
  lmer(Surface ~ (1 | genotype), data = seed_size_F9)
)$genotype %>%
  mutate(
    dummy_column = 0, # Include a column of zeroes to tell GEMMA to fit an intercept
    genotype   = row.names(.)
  ) %>%
  select(dummy_column, genotype, `(Intercept)`)
# Parents ar grouped into cohorts by year. INclude as a fixed effect
seed_size_parents_blups <- ranef(
  lmer(Surface ~ Year + (1 | genotype), data = seed_size_parents)
)$genotype %>%
  mutate(
    dummy_column = 0, # Include a column of zeroes to tell GEMMA to fit an intercept
    genotype   = row.names(.)
  ) %>%
  select(dummy_column, genotype, `(Intercept)`)

# Write files to disk
# For F9s, do this separately for replicates 1 and 2.
seedsize_F9_blups %>%
  filter(grepl("rep1", genotype)) %>%
  write_tsv(
    "03_processing/06_process_phenotypes/output/seed_size_blups_F9_rep1.tsv",
    col_names= FALSE
    )
seedsize_F9_blups %>%
  filter(grepl("rep2", genotype)) %>%
  write_tsv(
    "03_processing/06_process_phenotypes/output/seed_size_blups_F9_rep2.tsv",
    col_names= FALSE
  )
# Parents
seed_size_parents_blups %>%
  write_tsv(
    "03_processing/06_process_phenotypes/output/seed_size_blups_parents.tsv",
    col_names = FALSE
  )
