#' Rename BAM files
#'
#' List files with hard-to-read names including NGS facility info and barcodes,
#' create a tibble that lines them up with biological sample names, then rename
#' files accordingly.
#'
#' Input: Directory of BAM files
#' Output: The same BAM files, with new names.
#'
#' Tom Ellis, adapting code by Pieter Clauw, 25th November 2025

library(tidyverse)

# Location of aligned bams
bam_dir <- '/scratch-cbe/users/thomas.ellis/crosses/04_aligned_bam/'
bam_files <- list.files(path = bam_dir, pattern = '*.sort.bam$')
# Tible giving directory, NGS sample ID, and barcode sequence.
bam_info <- tibble(
  'directory' = bam_dir,
  'filename' = bam_files
) %>%
  separate(
    filename, sep = '_', into = c('NGS_sample_ID', 'barcode', NA, NA), remove = F
  )

# Details of sequencing plates
plate_info <- read_csv("01_data/02_F8_unaligned_bams/sequencing_plates_pieter.csv")

# Line up sample IDs with BAM-file names
plate_info <- plate_info %>%
  rename_with(~ gsub(' ', '_', .x)) %>% # Change spaces to underscores in column names
  mutate(
    name = gsub(' ', '_', name), # Change spaces to underscores in the variable `name`
    barcode = paste0(barcode2, barcode1),
    NGS_sample_ID = as.character(NGS_sample_ID)
    ) %>%
  left_join(bam_info, by = c('NGS_sample_ID', 'barcode'))


#' List entries with NA in the file path
#' Most are rows 7 to 12 which were used for other projects
#' The remaining three are real samples for which there is no library available.
plate_info %>%
  filter(is.na(directory))
# Filter them out.
plate_info <- plate_info %>%
  filter( !is.na(directory) )

# Rename BAM files to use the sample name.
file.rename(
  from = paste0(plate_info$directory, plate_info$filename),
  to   = paste0(plate_info$directory, plate_info$name, ".bam")
)
file.rename(
  from = paste0(plate_info$directory, plate_info$filename, '.bai'),
  to   = paste0(plate_info$directory, plate_info$name, ".bam.bai")
)
