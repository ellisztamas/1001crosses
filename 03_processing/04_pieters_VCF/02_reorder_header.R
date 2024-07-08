#' Create the files needed to correct the sample headers in Pieter's VCF file.
#'
#' The headers in Pieter's bam file are full paths. I want to replace these with
#' sample names, and also exclude samples that Pieter was not able to validate.
#' This creates text files to do this, and the actual heavy lifting happens in
#' 01_reheader_vcf.sh
#'
#' Inputs:
#'    plate_info.csv: A field created by Pieter summarising sequencing information
#'        Importantly, the last column 'geno_select' contains his assessment of
#'        whether each sample should be included based on his SNPmatch results.
#'    vcf_header_to_change.txt: list of sample headers in Pieter's VCF file that
#'        include full paths to files.
#' Outputs:
#'    new_vcf_header.txt: list of new header names for each sample in vcf_header_to_change.txt
#'    samples_to_keep.txt: list of samples that could be validated
#'    parents_to_keep.txt: list of parents

library(tidyverse)

# === Inputs ===

# Pieter's file giving samples including how he checked them. See:
# crosses_continued/004.F8/002.mapping/002.scripts/001.gather_pheno_geno.Rmd
plate_info <- read_csv(
  "../crosses_continued/004.F8/002.mapping/001.data/plate_info.csv",
  show_col_types = FALSE
  )

outdir <- '/scratch-cbe/users/thomas.ellis/crosses/04_pieters_VCF/01_reheader_VCF'
dir.create(outdir, showWarnings = FALSE)

# Sample names from Pieters VCF files to be changed
old_header <- read_csv(
  paste0(outdir, "/vcf_header_to_change.txt"),
  col_names = "filename",
  col_types = 'c'
)

# === Script ===

# Reorder sample names match the old VCF header
new_header <- old_header %>%
  mutate(
    filename = basename(filename)
  ) %>%
  left_join(plate_info, by='filename') %>%
  mutate(
    Sample_name = ifelse(genoSelection == "yes", Sample_name, paste0("xx_", Sample_name)),
    Sample_name = str_replace(Sample_name, "_F8_", "_")
    ) %>%
  select(Sample_name)

new_header %>%
  write_delim(
    paste0(outdir, "/new_vcf_header.txt"), col_names = FALSE, delim = " "
    )

new_header %>%
  filter(
    !grepl("xx_", Sample_name),
    !grepl("Aa1", Sample_name)
    ) %>%
  write_delim(
    paste0(outdir, "/samples_to_keep.txt"), col_names = FALSE, delim = " "
  )

# Text file giving parents of any F8 in either replicate cohort.
new_header %>%
  filter(
    !grepl("xx_", Sample_name),
    !grepl("Aa1", Sample_name)
  ) %>%
  mutate(
    Sample_name = str_replace(Sample_name, "_rep[1,2]", "")
  ) %>%
  separate(Sample_name, into =c("p1", 'p2'), sep="x") %>%
  pivot_longer(p1:p2) %>%
  select(value) %>%
  distinct() %>%
  arrange(value) %>%
  write_csv(
    paste0(outdir, "/parents_to_keep.txt"), col_names = FALSE
  )

