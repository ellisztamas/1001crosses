
library(tidyverse)

snpmatch_files <- Sys.glob("03_processing/pieters_sample_sheet/output/snpmatch/*.csv")
snpmatch_files <- snpmatch_files[grep("_rep[12].csv", snpmatch_files)]

sm_results <- lapply(snpmatch_files, read_csv, show_col_types=FALSE) %>%
  do.call(what = 'rbind')

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

write_csv(
  sm_results,
  "03_processing/pieters_sample_sheet/output/snpmatch/SNPmatch_summary.csv"
  )

sm_results %>%
  filter(CorrectParents == "no") %>%
  select(FOL, ExpectedParent1, ExpectedParent2, ObservedParent1, ObservedParent2)





