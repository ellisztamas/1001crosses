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

# Data frame giving original and corrected names
source("03_processing/06_process_phenotypes/01_update_names.R")

# List of parental and F8 lines that have been genotyped and validated
parental_names <- read_csv(
  "03_processing/09_impute_haplotypes/output/parental_line_names.txt",
  col_names = c("genotype"), col_types = 'c'
)
progeny_names <- read_csv(
  "03_processing/09_impute_haplotypes/output/F8_imputed_line_names.txt",
  col_names = c("genotype"), col_types = 'c'
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

# Swap genotype names for corrected names where necessary
flowering_time <-flowering_time %>%
  left_join(naming_file, by=c("genotype" = "original_name")) %>%
  mutate(
    corrected_name = ifelse(generation == "parent", genotype, corrected_name)
  ) %>%
  select(-genotype) %>%
  dplyr::rename(genotype = corrected_name)

# BLUPs for flowering time, accounting for tray nested within replicate, and
# cohort (parent, F8 cohort 1, F8 cohort 2).
mod_ft <- glmer(
  days_to_flower ~ (1 | genotype) + (1| replicate/tray) + (1| cohort),
  data=flowering_time, family = 'poisson'
)

# Residuals are somewhat overdispersed
# hist(fitted(mod_ft))
# hist(resid(mod_ft))
# plot(fitted(mod_ft), resid(mod_ft))

ft_blups <- ranef(mod_ft)$genotype %>%
  mutate(
    dummy_column = 0, # Include a column of zeroes to tell GEMMA to fit an intercept
    genotype   = row.names(.)
  ) %>%
  dplyr::rename(
    blup = `(Intercept)`
  )


# Data file for parental lines.
# 15 lines have no phenotype. They either didn't germinate or didn't flower
parental_names %>%
  left_join(ft_blups) %>% # Left join preserves the order in line names
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
  left_join(ft_blups) %>%
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
