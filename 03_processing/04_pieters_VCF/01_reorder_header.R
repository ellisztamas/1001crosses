#'
library(tidyverse)

# Pieter's file giving samples including how he checked them. See:
# crosses_continued/004.F8/002.mapping/002.scripts/001.gather_pheno_geno.Rmd
plate_info <- read_csv(
  "../crosses_continued/004.F8/002.mapping/001.data/plate_info.csv",
  show_col_types = FALSE
  )

# Sample names from Pieters VCF files to be changed
old_header <- read_csv(
  "03_processing/pieters_VCF/vcf_header_to_change.txt",
  col_names = "filename",
  col_types = 'c'
)

# Reorder sample names match the old VCF header
new_header <- old_header %>%
  mutate(
    filename = basename(filename)
  ) %>%
  left_join(plate_info, by='filename') %>%
  mutate(
    Sample_name = ifelse(genoSelection == "yes", Sample_name, paste0("xx_", Sample_name))
    ) %>%
  select(Sample_name)

new_header %>%
  write_delim(
    "03_processing/pieters_VCF/new_vcf_header.txt", col_names = FALSE, delim = " "
    )

new_header %>%
  filter(!grepl("xx_", Sample_name)) %>%
  write_delim(
    "03_processing/pieters_VCF/samples_to_keep.txt", col_names = FALSE, delim = " "
  )