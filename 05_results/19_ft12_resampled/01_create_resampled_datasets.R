#' Create new F9 GWAS populations
#'
#' There are two independent replicate populations of F9s.
#' Most lines have a matching pair in each cohort.
#'
#' This script creates new datasets of roughly the same size, ensuring that
#' at most line one of each cross is included.
#'
#' Tom Ellis, 23rd June 2025


set.seed(213)

library(tidyverse)


# Input -------------------------------------------------------------------

ft12 <- read_delim(
  "03_processing/06_process_phenotypes/output/flowering_time_blups_F9_combined.tsv",
  delim="\t",
  col_names = c("FID", "ID", "phenotype"),
  col_types = 'ccd'
)

# Output ------------------------------------------------------------------

outdir <- "/scratch-cbe/users/thomas.ellis/crosses/19_ft12_resampled/"


# Main --------------------------------------------------------------------

# Vector of unique cross names without cohort, e.g. "9391x9353"
cross_list <- ft12 %>%
  separate(ID, into = c("cross", "cohort"), sep = "_", remove = FALSE) %>%
  pull(cross) %>%
  unique()

# Create 200 new datasets
for(i in 1:200){
  # Vector of zero and ones.
  cohort_01 <- rbinom(length(cross_list), size=1, p=0.5)
  # New line names, with either '_rep1' or "_rep2"
  new_cohort_1 <- paste0(cross_list, "_rep", cohort_01 + 1)
  # Create a new dataset
  new_dataset <- ft12 %>%
    filter(
      ID %in% new_cohort_1, # Keep only those lines in new_cohort_1
      !is.na(phenotype)
    ) %>%
    mutate(
      phenotype = phenotype - mean(phenotype),
      phenotype = phenotype / sd(phenotype)
    )
  cat('Dataset', i, 'created with', nrow(new_dataset), "lines.\n")

  # Double check that there are no duplicates of any line
  line_counts <- new_dataset %>%
    separate(ID, into = c("cross", "cohort"), sep = "_", remove = FALSE) %>%
    pull(cross) %>%
    table()
  if(any(line_counts != 1)) {
    cat("One or more lines are replicated more than once at iteration", i, "\n")
  }

  # Save to disk
  dir.create(
    paste0(outdir, "ft12_resample_", i),
    showWarnings = FALSE, recursive=TRUE
    )
  new_dataset %>%
    write_tsv(
      paste0(outdir, "ft12_resample_", i, "/phenotypes.tsv"),
      col_names = FALSE
    )
}
