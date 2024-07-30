#' Prepare raw rosette size data for further analysis.
#'
#' Estimate genetic values for each line from the phenotype data by fitting
#' BLUPs.
#'
#' Note that accession 1137 does not appear to have been sequenced, so this is
#' removed for now.
#'
#' Inputs:
#'    Raw data file for rosette size including replicate lines 1 and 2 and parents
#'    from the growth-room phenotyping experiment.
#'
#' Outputs:
#'    Tab-separated text files formatted for GEMMA.
#'    The first column is a vector of 1s to tell GEMMA to fit an intercept.
#'    Subsequent columns show genotype and BLUP estimates of
#'    genetic values of rosette size. There are separate file for replicate lines
#'    1 and 2, plus the parents

library(tidyverse)
library(lme4)

source("03_processing/04_pieters_VCF/04_align_original_corrected_sample_names.R")

# Import data on rosette size for each individual
raw_rosette_size <- read_csv("01_data/05_phenotype_expt/rosette_size_parents_F9s.csv") %>%
  filter(
    cohort != 1, # There are no data on cohort 1
    size1 > 0 # Remove plants that never germinated.
    ) %>%
  # Rosette size is the mean of three measurements
  mutate(
    genotype = str_replace(genotype, " ", "_"),
    rosette = (size1 + size2 + size3) / 3
  )

# Correct sample names and remove dubious genotypes.
raw_rosette_size <- raw_rosette_size %>%
  left_join(sample_name_corrections, by=c("genotype" = "Sample_name_original")) %>%
  # If the sample is a parent, add genotpe to Sample_name_corrected
  mutate(
    Sample_name_corrected = ifelse(generation == "parent", genotype, Sample_name_corrected)
  ) %>%
  filter( ! is.na(Sample_name_corrected)) %>% # Remove entries that are still NA
  rename(
    genotype_original = genotype,
    genotype = Sample_name_corrected
  )

# Add a column 'cohort' indicating whether a sample is a parent, or from cross
# replicate1 or 2
raw_rosette_size <- raw_rosette_size %>%
  mutate(
    cohort = generation,
    cohort = ifelse(grepl("rep1", genotype), "rep1", generation),
    cohort = ifelse(grepl("rep2", genotype), "rep2", generation)
  )

# Remove 5835, because Pieter thinks this was confused with 6180
# Also 1435, which Pieter think was confused with 7383
raw_rosette_size <- raw_rosette_size %>%
  filter(
    genotype != "5835",
    genotype != "1435"
    )

# Fit BLUPs for rosette size
# Notice the fixed effect of 'generation' which allows for different means for
rosette_blups <- ranef(
  lmer( log(rosette) ~ cohort + (1 | tray) + (1 | genotype), data = raw_rosette_size)
)$genotype %>%
  mutate(
    dummy_column = 0, # Include a column of zeroes to tell GEMMA to fit an intercept
    genotype   = row.names(.)
  ) %>%
  select(dummy_column, genotype, `(Intercept)`)

rosette_blups %>%
  filter(
    grepl("rep1", genotype)
  ) %>%
  write_tsv(
    "03_processing/06_process_phenotypes/output/rosette_size_blups_F9_rep1.tsv",
    col_names = FALSE
    )

rosette_blups %>%
  filter(
    grepl("rep2", genotype)
  ) %>%
  write_tsv(
    "03_processing/06_process_phenotypes/output/rosette_size_blups_F9_rep2.tsv",
    col_names = FALSE
  )

rosette_blups %>%
  filter(
    ! grepl("rep", genotype),
    genotype != "1137"
  ) %>%
  write_tsv(
    "03_processing/06_process_phenotypes/output/rosette_size_blups_parents.tsv",
    col_names = FALSE
  )