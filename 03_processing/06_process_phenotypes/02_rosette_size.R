#' Prepare raw rosette size data for further analysis.
#'
#' Estimate genetic values for each line from the phenotype data by fitting
#' BLUPs.
#'
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
#'    1 and 2, plus the parents.
#'    Output files are in the same order as the genotype files, with missing
#'    phenotype values where necessary.

library(tidyverse)
library(lme4)

# === Inputs === #

#' Raw data file for rosette size including replicate lines 1 and 2 and parents
#' from the growth-room phenotyping experiment.
raw_data_path <- "01_data/05_phenotype_expt/rosette_size_parents_F9s.csv"

# Data frame giving original and corrected names
source("03_processing/06_process_phenotypes/01_update_names.R")
# List of parental and F8 lines that have been genotyped and validated
parental_names <- read_csv(
  "03_processing/05_imputation/output/parental_line_names.txt",
  col_names = c("genotype"), col_types = 'c'
)
progeny_names <- read_csv(
  "03_processing/05_imputation/output/F8_phased_imputed_line_names.txt",
  col_names = c("genotype"), col_types = 'c'
)



# === Output === #

# Directory for saving the output files (which will be the input for GEMMA)
outdir <- "03_processing/06_process_phenotypes/output"
dir.create(outdir, showWarnings = FALSE)

# Text files for three cohorts that PLINK can read
rep1    <- paste0(outdir, "/rosette_size_blups_rep1.tsv")
rep2    <- paste0(outdir, "/rosette_size_blups_rep2.tsv")
combined<- paste0(outdir, "/rosette_size_blups_F9_combined.tsv")
parents <- paste0(outdir, "/rosette_size_blups_parents.tsv")



# === Main === #

# Import data on rosette size for each individual
raw_rosette_size <- read_csv(raw_data_path) %>%
  rename(
    replicate = cohort
  ) %>%
  filter(
    replicate != 1, # There are no data on cohort 1
    size1 > 0 # Remove plants that never germinated.
    ) %>%
  mutate(
    genotype = str_replace(genotype, " ", "_")
  )
# Swap genotype names for corrected names where necessary
raw_rosette_size %>%
  left_join(naming_file, by =c("genotype" = "original_name")) %>%
  mutate(
    corrected_name = ifelse(generation == "parent", genotype, corrected_name)
  ) %>%
  select(-genotype) %>%
  rename(genotype = corrected_name)
# Add a column 'cohort' indicating whether a sample is a parent, or from cross
# replicate1 or 2
raw_rosette_size <- raw_rosette_size %>%
  mutate(
    cohort = case_when(
      generation == "parent" ~ "parent",
      grepl("rep1", genotype) ~ "rep1",
      grepl("rep2", genotype) ~ "rep2"
    )
  )
# Pivot longer so that all three size estimates for each plant are on separate rows.
raw_rosette_size <- raw_rosette_size %>%
  pivot_longer(size1:size3, names_to = "measurement", values_to = "rosette")

# Square-root transforming the data gives surprisingly nicely normal data
# hist(raw_rosette_size$rosette)
# hist(sqrt(raw_rosette_size$rosette))
# hist(log(raw_rosette_size$rosette))

# Fit BLUPs for rosette size
rosette_model <- lmer(
    rosette ~ (1|cohort) + (1 | tray) + (1 | genotype/id) + (1| replicate) ,
    data = raw_rosette_size
    )

# The residuals are overdispersed!
# This doesn't change if I use log or sqrt transforms (sqrt helps, but reduces interpretebility)
# The issue goes away if I do:
# sqrt(rosette) ~ (1|cohort) + (1 | tray) + (1 | genotype)
# However,  plant ID seems to explain a lot of variation
# hist(fitted(rosette_model))
# hist(resid(rosette_model))
# plot(fitted(rosette_model), resid(rosette_model))

rosette_blups <- ranef(rosette_model)$genotype %>%
  mutate(
    dummy_column = 0, # Include a column of zeroes to tell GEMMA to fit an intercept
    genotype   = row.names(.)
  ) %>%
  select(dummy_column, genotype, `(Intercept)`) %>%
  rename(
    blup = `(Intercept)`
  )

# Data file for parental lines.
# 15 lines have no phenotype. They either didn't germinate or didn't flower
parental_names %>%
  left_join(rosette_blups) %>% # Left join preserves the order in line names
  select(dummy_column, genotype, blup) %>%
  mutate(
    dummy_column = 0
  ) %>%
  filter(
    !is.na(blup)
  ) %>%
  write_tsv(parents, col_names = FALSE)

# BLUPs for progeny, in the order they appear in the genotype file.
progeny_blups <- progeny_names %>%
  left_join(rosette_blups) %>%
  select(dummy_column, genotype, blup) %>%
  mutate(
    dummy_column = 0
    ) %>%
  filter(
    !is.na(blup)
  )
# Data file for all progeny lines combined
progeny_blups %>%
  write_tsv(combined, col_names = FALSE)
# Data file for all cohort 1
progeny_blups %>%
  filter(
    grepl("rep1", genotype)
  ) %>%
  write_tsv(rep1, col_names = FALSE)
# Data file for all cohort 2
progeny_blups %>%
  filter(
    grepl("rep2", genotype)
  ) %>%
  write_tsv(rep2, col_names = FALSE)

