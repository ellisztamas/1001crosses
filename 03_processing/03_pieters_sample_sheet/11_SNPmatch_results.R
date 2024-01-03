
library(tidyverse)

# Import the output of 02_library/CSV_from_SNPmatch.py for each sample, and merge
snpmatch_files <- Sys.glob("03_processing/03_pieters_sample_sheet/output/snpmatch/*.csv")
snpmatch_files <- snpmatch_files[grep("_rep[12].csv", snpmatch_files)] # ignore ancestry files.
sm_results <- lapply(snpmatch_files, read_csv, show_col_types=FALSE) %>%
  do.call(what = 'rbind')
# Create a column stating whether the observed parents match the expected parents
sm_results <- sm_results %>%
  separate(
    ExpectedParents,
    sep = 'x', into = c('ExpectedParent1', 'ExpectedParent2')
    ) %>%
  mutate(
    CorrectParents = ifelse(
      (ExpectedParent1 == ObservedParent1 | ExpectedParent1 == ObservedParent2) &
        (ExpectedParent2 == ObservedParent1 | ExpectedParent2 == ObservedParent2),
      'yes', 'no')
    )

# write_csv(
#   sm_results,
#   "03_processing/03_pieters_sample_sheet/output/snpmatch/SNPmatch_summary.csv"
#   )

sm_results %>%
  filter(CorrectParents == "no") %>%
  select(FOL, ExpectedParent1, ExpectedParent2, ObservedParent1, ObservedParent2)

# about 2/3 crosses have two replicates
sm_results %>%
  filter(CorrectParents == "yes") %>%
  separate(FOL, into=c("cross", 'generation', 'rep')) %>%
  group_by(cross) %>%
  summarise(
    n = n()
  ) %>%
  pull(n) %>%  table

# 167 unique crosses, of which two are reciprocal
sm_results %>%
  # filter(CorrectParents == "yes") %>%
  separate(FOL, into=c("cross", 'generation', 'rep')) %>%
  mutate(reciprocal = paste0(ExpectedParent2, "x", ExpectedParent1)) %>%
  select(cross, reciprocal) %>%
  unique() %>%
  filter(reciprocal %in% cross)



