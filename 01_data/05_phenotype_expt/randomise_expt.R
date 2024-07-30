#' Randomise F9 and parental plants to quantify flowering time.
#'
#' The goal is to grow three replicates each of 427 F9 families with the 219
#' parents in three fully randomised cohorts. Each tray can hold 48 pots, meaning
#' we need 13 full trays and half of a 14th to hold a single replicate cohort.
#'
#' This script takes a list of F9 families for which seed is available (created
#' by Tal Dahan for her seed size project) and extracts the IDs of the parents of
#' those families.
#'
#' Input:
#'    CSV of available seed for F9 families created
#' Output:
#'    CSV giving design of the experiment. Columns show cohort (1-3), tray (1-42),
#'    row within each tray (A-H), column within each tray (1-6) and genotype for
#'    each pot.
#'
#' Tom Ellis, 5th February 2024

library(tidyverse)

set.seed(646)

# Vector giving IDs of all the F8s
F8_names <- read_csv(
  "01_data/05_phenotype_expt/F9_seed_availability.csv",
  col_types = 'ccc') %>%
  pull(id)

# Extract a list of names for the parents from that.
parent_names <- F8_names %>%
  str_extract("[0-9]{3,4}x[0-9]{3,4}") %>%
  str_split("x") %>%
  do.call(what = 'c') %>%
  unique()

n_lines <- length(F8_names) + length(parent_names)


# Dataframes giving cohort, tray, row and column for each of 14 complete trays
# of 48 plants in three cohorts.
ft10_design <- vector('list', 3)
for(cohort in 1:3){
  this_cohort <- vector('list', 14)
  for(tray in 1:14){
    this_cohort[[tray]] <-  tibble(
      cohort = cohort,
      tray = tray + (cohort-1)*14,
      row = rep(LETTERS[1:8], 6),
      col = rep(1:6, each=8)
    )
  }
  ft10_design[[cohort]] <- do.call('rbind', this_cohort)[1:n_lines,]
  # Add genotypes
  ft10_design[[cohort]]$genotype = sample(
    c(parent_names, F8_names), size = n_lines, replace = FALSE
  )
}
ft10_design <- do.call('rbind', ft10_design) %>%
  mutate(
    generation = ifelse(grepl(x=genotype, pattern="x"), "F8", "parent")
  ) %>%
  arrange(generation, genotype) %>%
  mutate(
    id  = 1:nrow(.)
  ) %>%
  arrange(cohort, tray, row, col) %>%
  select (id, cohort, tray, row, col, genotype, generation)

# As a sanity check, ensure all genotypes are replicated exactly 3 times.
all(
  table(ft10_design$genotype) == 3
  )

write_csv(ft10_design, file = "01_data/05_phenotype_expt/randomise_expt.csv")

ft10_design %>%
  mutate(
    pos = paste(row, col, sep = ""),
    pos = paste(cohort,tray, pos, sep = ".")
    ) %>%
  arrange(generation, genotype) %>%
  select(id, genotype, pos) %>%
  write_tsv(
    "01_data/05_phenotype_expt/randomise_expt.txt",
    col_names = FALSE
    )
