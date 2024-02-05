
library(tidyverse)

# Import the output of 02_library/CSV_from_SNPmatch.py for each sample, and merge
snpmatch_files <- Sys.glob("03_processing/03_pieters_sample_sheet/output/snpmatch/*.csv")
snpmatch_files <- snpmatch_files[grep("_rep[12].csv", snpmatch_files)] # ignore ancestry files.
sm_results <- lapply(snpmatch_files, read_csv, col_types = cols(.default = "c")) %>%
  do.call(what = 'rbind')
# Create a column stating whether the observed parents match the expected parents
sm_results <- sm_results %>%
  select(FOL, ExpectedParents, ObservedParent1, ObservedParent2, HomozygousWindow) %>%
  separate(FOL, into=c("cross", 'generation', 'rep')) %>%
  group_by(cross) %>%
  separate(
    ExpectedParents,
    sep = 'x', into = c('ExpectedParent1', 'ExpectedParent2')
    ) %>%
  mutate(
    top_match = ifelse(
      (ExpectedParent1 == ObservedParent1) & (ExpectedParent2 == ObservedParent2) |
        (ExpectedParent1 == ObservedParent2) & (ExpectedParent2 == ObservedParent1),
      "yes", "no"
    )
  )


#' List of vectors, each containing the IDs of compatible parents for each sample
list_of_matches <- sm_results$HomozygousWindow %>% str_split(",")
# Add columns indicating whether each expected parent is among the compatible candidates
sm_results$found_parent1 <- NA
sm_results$found_parent2 <- NA
for(i in 1:nrow(sm_results)){
  sm_results$found_parent1[i] <- sm_results$ExpectedParent1[i] %in% list_of_matches[[i]]
  sm_results$found_parent2[i] <- sm_results$ExpectedParent2[i] %in% list_of_matches[[i]]
}
sm_results$found_both_parents <- sm_results$found_parent1 & sm_results$found_parent2

# write_csv(
#   sm_results,
#   "03_processing/03_pieters_sample_sheet/output/snpmatch/SNPmatch_summary.csv"
#   )

sm_results %>%
  filter(CorrectParents == "no") %>%
  select(FOL, ExpectedParent1, ExpectedParent2, ObservedParent1, ObservedParent2) %>%
  View()

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



