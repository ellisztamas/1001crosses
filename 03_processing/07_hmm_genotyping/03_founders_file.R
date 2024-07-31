library(tidyverse)

read_csv(
  "03_processing/04_pieters_VCF/output/samples_to_keep.txt",
  col_names = "offspring") %>%
  # filter(
  #   grepl("rep1", offspring)
  # ) %>%
  mutate(
    parents = str_replace(offspring, "_rep.", "")
  ) %>%
  separate(parents, into = c("parent1", 'parent2'), sep="x") %>%
  write_delim(
    "03_processing/07_hmm_genotyping/output/founders.tsv", col_names = FALSE
  )
