#' Create a list
#'
#' I received seed bags back from Labdeers in a specific order.
#' This script aligns the position labels from those bags with their legacy
#' genotypes from before validation, and also assign a new genotype label.
#' The new genotype label is in the format ISC_345.1 where
#'     - ISC is intra-Swedish crosses
#'     - 145 is a unique integer
#'     - .1 or .2 indicates whether the cross was originally part of cross cohort
#'         1 or cohort 2.
#'
#' This script only applies to lines that generated *at least one* seed bag in
#' the phenotyping experiment.
#'
#' Tom Ellis, 3rd September 2025

library("tidyverse")

catalogue <- read_csv(
  "01_data/10_grow_F10s/F10_seed_size_catalogue_labdeers.csv"
  )

catalogue <- catalogue %>%
  mutate(
    label = ifelse(is.na(label_labdeers), label, label_labdeers)
  )

exp_design <- read_tsv(
  "01_data/05_phenotype_expt/randomise_expt.txt",
  col_names = c('i', 'genotype', 'label')
)

exp_design <- exp_design %>%
  filter(grepl(" rep", genotype))

new_names <- exp_design %>%
  select(genotype) %>%
  distinct() %>%
  mutate(
    id = paste0(
      "ISC_",
      str_pad(1:nrow(.), 3, pad="0"),
      ".",
      str_sub(genotype, -1,-1)
    )
  )


catalogue %>%
  full_join(exp_design, by = 'label') %>%
  select(label, genotype) %>%
  full_join(
    new_names,
    by = 'genotype') %>%
  filter(!is.na(genotype)) %>%
  mutate(
    pool = rep(1:4, nrow(.)/4),
    i = 1:1284
  ) %>%
  write_tsv(
    "01_data/11_field_experiment/seed_tube_checklist.tsv"
    )



# check <- read_tsv(
#   "01_data/10_grow_F10s/seed_bag_checklist.tsv"
# )
# check %>%
#   group_by(genotype) %>%
#   filter(n() > 1)
# check %>%
#   group_by(id) %>%
#   filter(n() > 1)


