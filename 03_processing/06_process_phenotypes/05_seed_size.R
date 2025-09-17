library('lme4')
library(tidyverse)

# === Inputs === #

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

# File giving genotype and position in the phenotyping experiment
exp_design <- read_tsv(
  "01_data/05_phenotype_expt/randomise_expt.txt",
  col_names = c('i', 'genotype', 'label'),
  show_col_types = FALSE
)

# Raw data on seed size from Labdeers, with a row for each seed
raw_seed_sizes <- read_delim(
  "01_data/05_phenotype_expt/seed_sizes_per_seed.txt",
  show_col_types = FALSE
)


# === Output ====

# Directory for saving the output files (which will be the input for GEMMA)
outdir <- "03_processing/06_process_phenotypes/output"
dir.create(outdir, showWarnings = FALSE)

# Text files for three cohorts that PLINK can read
rep1    <- paste0(outdir, "/seed_size_blups_rep1.tsv")
rep2    <- paste0(outdir, "/seed_size_blups_rep2.tsv")
combined<- paste0(outdir, "/seed_size_blups_F9_combined.tsv")
parents <- paste0(outdir, "/seed_size_blups_parents.tsv")



# === Main ===

# Return row ID, genotype, label, replicate, tray and position
exp_design <- exp_design %>%
  separate(
    col = 'label',
    into = c('rep', 'tray', 'pos'),
    remove = FALSE
  ) %>%
  # Correct line names
  mutate(
    genotype = gsub(" rep", "_rep", genotype)
  )

# Simplify to label, intrument, pixel count, sirface area and SSE, then
# merge with exp_design
raw_seed_sizes <- raw_seed_sizes %>%
  dplyr::select(`Pick description`, Instrument, `Pixel count`, Surface, SSE) %>%
  dplyr::rename(
    label=`Pick description`,
    instrument = Instrument,
    pixel_count = `Pixel count`,
    surface_area = Surface
  ) %>%
  # Merge with exp_design to get genotypes and positions.
  left_join(exp_design, by = 'label') %>%
  # Correct line names
  left_join(naming_file, by=c("genotype" = "original_name")) %>%
  mutate(
    corrected_name = ifelse(grepl("x", genotype), corrected_name, genotype),
    cohort = case_when(
      !grepl("rep", genotype) ~ "parent",
      grepl("rep1", genotype) ~ "rep1",
      grepl("rep2", genotype) ~ "rep2"
    )
  ) %>%
  dplyr::select(-genotype) %>%
  dplyr::rename(genotype = corrected_name)

# Line means look mostly normal, with some right skew.
# raw_seed_sizes %>%
#   group_by(label) %>%
#   summarise(
#     surface_area = mean(surface_area)
#   ) %>%
#   ggplot(aes(x=surface_area)) +
#   geom_histogram()

# Variance decomposition.
# Label (i.e. plant ID) is required, because there are many seeds per plant
model_seed_size <- lmer(
  surface_area ~ (1 | rep/tray) + (1 | genotype/label) + (1|cohort),
  data = raw_seed_sizes
  )

# As ever, residuals look overdispersed.
# hist(fitted(model_seed_size))
# hist(resid(model_seed_size))
# plot(fitted(model_seed_size), resid(model_seed_size))

seed_size_blups <- ranef(model_seed_size)$genotype %>%
  mutate(
    dummy_column = 0, # Include a column of zeroes to tell GEMMA to fit an intercept
    genotype   = row.names(.)
  ) %>%
  dplyr::rename(
    blup = `(Intercept)`
  )

# Residuals look fine
# Note that using log seed size introduces a pattern; raw values are better
# plot(fitted(model_seed_size), resid(model_seed_size))

# There are a few outliers.
# hist(ranef(model_seed_size)$genotype[,1])
# ranef(model_seed_size)$genotype %>%
#   filter(`(Intercept)` > 0.03)


# Data file for parental lines.
# 11 lines have no phenotype. They either didn't germinate or didn't flower
parental_names %>%
  left_join(seed_size_blups) %>% # Left join preserves the order in line names
  select(dummy_column, genotype, blup) %>%
  mutate(
    dummy_column = 0
  ) %>%
  filter(
    !is.na(blup)
  ) %>%
  write_tsv(parents, col_names = FALSE)


# BLUPs for progeny, in the order they appear in the genotype file.
# 12 lines have no data.
progeny_blups <- progeny_names %>%
  left_join(seed_size_blups) %>%
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

