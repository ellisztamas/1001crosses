#' Prepare line-mean flowering time data for GEMMA
#'
#' Calculates BLUPs for flowering time accounting for the experimental structure,
#' splits these into rep1, rep2 and parents. The output is saved in a format
#' PLINK can use, including a column of zeroes for the intercept.
#'
#'
#' Tom Ellis, 17th January 2025

library('tidyverse')
library('lme4')

# === Inputs === #

# Import tidied flowering time data
source("01_data/05_phenotype_expt/format_flowering_time_data.R")

# Genotyping validation results
ibdpainting <- read_csv(
  "03_processing/03_validate_genotypes/output/ibdpainting_results.csv"
)


# === Output === #

# Directory for saving the output files (which will be the input for GEMMA)
outdir <- "03_processing/06_process_phenotypes/output"
dir.create(outdir, showWarnings = FALSE)

# Text files for three cohorts that PLINK can read
rep1    <- paste0(outdir, "/flowering_time_blups_rep1.tsv")
rep2    <- paste0(outdir, "/flowering_time_blups_rep2.tsv")
combined<- paste0(outdir, "/flowering_time_blups_F9_combined.tsv")
parents <- paste0(outdir, "/flowering_time_blups_parents.tsv")

# === Main === #

# BLUPs for flowering time.
mod_ft <- glmer(
  days_to_flower ~ (1 | genotype) + (1| cohort/tray) + (1| genotype:cohort),
  data=flowering_time, family = 'poisson'
)
ft_blups <- ranef(mod_ft)$genotype %>%
  mutate(
    dummy_column = 0, # Include a column of zeroes to tell GEMMA to fit an intercept
    genotype   = row.names(.)
  ) %>%
  select(dummy_column, genotype, `(Intercept)`)


# Filter out lines that don't have good genotype data
lines_to_exclude <- ibdpainting %>%
  mutate(
    name = gsub(" F8 ", "_", name)
    ) %>%
  filter(
    grepl("no data|Self", diagnosis)
  ) %>%
  pull(name)
# Add two parental lines still be to sequenced
lines_to_exclude <- c(lines_to_exclude, "1137", "1074")
# Exclude them from the phenotype object
ft_blups <- ft_blups %>%
  filter(
    !genotype %in% lines_to_exclude
  )




# Create data files
# Values for F9 rep 1
ft_blups %>%
  filter(
    grepl("rep1", genotype)
    ) %>%
  write_tsv(rep1, col_names = FALSE)

# Values for F9 rep2
ft_blups %>%
  filter(
    grepl("rep2", genotype)
  ) %>%
  write_tsv(rep2, col_names = FALSE)

# Values for the two F9 cohorts combined
ft_blups %>%
  filter(
    grepl("rep[12]", genotype)
    ) %>%
  write_tsv(combined, col_names = FALSE)


# Values for the parents
ft_blups %>%
  filter(
    !grepl("rep[12]", genotype)
    ) %>%
  write_tsv(parents, col_names = FALSE)
