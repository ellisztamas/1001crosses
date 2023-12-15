#' Prepare flowering-time data the parents for GWAS.
#'
#' Import data on number of days to flowering at 10°C for the parents, and the
#' phenotype file previously created for the F8s with confirmed genotypes.
#'
#' Input: Flowering-time data for the F8s and parents
#' Output: PLINK formatted phenotype and covariate files for GEMMA
#'
#' Tom Ellis, 12th December 2023.

parents <- read_csv(
  "01_data/05_phenotypes/flowering_time_10C_parents.csv",
  col_types="cd"
  )

# f8s <-
read_tsv(
  "05_results/01_ft10_offspring/input/ft10_rep1.tsv",
  col_names = c("null_column","sample", "FT10C"),
  col_types = 'icd'
  ) %>%
  select(sample, FT10C) %>%
  separate(sample, sep="_", into=c("cross", "generation", "rep")) %>%
  separate(cross, into=c("parent1", 'parent2'), sep="x") %>%
  pivot_longer(parent1:parent2) %>%
  pull(value) %>%
  unique()

