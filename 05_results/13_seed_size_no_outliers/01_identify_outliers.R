#' Identify and remove accessions with unusually large seeds.

library(tidyverse)


# === Input files === #

# Import original BLUP files
filenames <- list(
  parents = "03_processing/06_process_phenotypes/output/seed_size_blups_parents.tsv",
  rep1    = "03_processing/06_process_phenotypes/output/seed_size_blups_rep1.tsv",
  rep2    = "03_processing/06_process_phenotypes/output/seed_size_blups_rep2.tsv"
)
seed_size_blups <- lapply(
  filenames, read_delim,
  delim = "\t",
  col_names = c('intercept', 'genotype', 'phenotype'), col_types = 'icd'
)

# === Output === #

outdir <- "05_results/13_seed_size_no_outliers/output"
dir.create(outdir, showWarnings = FALSE)

# === Script === #

# There is not really much skew among the parents
seed_size_blups$parents %>%
  ggplot(aes(x = phenotype)) +
  geom_histogram()
# There are 3 lines that look outlier-y in cohort 1
seed_size_blups$rep1 %>%
  ggplot(aes(x = phenotype)) +
  geom_histogram()
# Rep2 has 2
seed_size_blups$rep2 %>%
  ggplot(aes(x = phenotype)) +
  geom_histogram()

# There is not much overlap between the biggest genotypes in the three cohorts
# Only 6114x9381 appears twice
seed_size_blups$parents %>%
  filter(
    phenotype > 0.025
  )
seed_size_blups$rep1 %>%
  filter(
    phenotype > 0.025
  )
seed_size_blups$rep2 %>%
  filter(
    phenotype > 0.025
  )


# Filter them out anyway and save for GWAS
for(cohort in names(seed_size_blups)){
  outfile = paste0(
    outdir, "/seed_size_", cohort, "_no_outliers.tsv"
  )
  seed_size_blups[[cohort]] %>%
    filter(
      phenotype < 0.025
    ) %>%
    write_delim(outfile, delim = "\t", col_names = FALSE)
}
