#' Identify and remove accessions with unusually large seeds.

library(tidyverse)

Sys.glob("03_processing/06_process_phenotypes/output/seed_size_blups_*.tsv")

# === Input files === #

# Import original BLUP files
filenames <- list(
  parents = "03_processing/06_process_phenotypes/output/seed_size_blups_parents.tsv",
  rep1    = "03_processing/06_process_phenotypes/output/seed_size_blups_F9_rep1.tsv",
  rep2    = "03_processing/06_process_phenotypes/output/seed_size_blups_F9_rep2.tsv"
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

# Parents and rep1 have three outliers with BLUP values > 0.03
seed_size_blups$parents %>%
  ggplot(aes(x = phenotype)) +
  geom_histogram()
seed_size_blups$rep1 %>%
  ggplot(aes(x = phenotype)) +
  geom_histogram()
# Rep2 has 1
seed_size_blups$rep2 %>%
  ggplot(aes(x = phenotype)) +
  geom_histogram()

# There is not much overlap between the biggest genotypes in the three cohorts
seed_size_blups$parents %>%
  filter(
    phenotype > 0.03
  )
seed_size_blups$rep1 %>%
  filter(
    phenotype > 0.03
  )
seed_size_blups$rep2 %>%
  filter(
    phenotype > 0.03
  )


# Filter them out anyway and save for GWAS
for(cohort in names(seed_size_blups)){
  outfile = paste0(
    outdir, "/seed_size_", cohort, "_no_outliers.tsv"
  )
  seed_size_blups[[cohort]] %>%
    filter(
      phenotype < 0.03
    ) %>%
    write_delim(outfile, delim = "\t", col_names = FALSE)
}
