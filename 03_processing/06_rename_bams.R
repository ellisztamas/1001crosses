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
library(googlesheets4)
gs4_deauth()

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

# Plates 1-4 are in one google sheet, and plate 5 in another.
# Note that Pieter Clauw found an error in the sheet for plate 5 (the plate was
# processed upside down), which he manually corrected.
sheets_list <- list(
  plates1_4 = "https://docs.google.com/spreadsheets/d/16hyDkyUjRDIQiLbczuV986bWgC3MkNkHrtqGkJuUZs4/edit#gid=1348764971",
  plate5    = "https://docs.google.com/spreadsheets/d/1ZnW0hLDiGwWzfNsCf-gubHBbM_T8H4F0svZGp47dvU8/edit#gid=0"
)
# Import the data sheets
plate_info <- vector('list', 5)
for(plate in 1:5){
  # Decide whether the plate is 1-4 or 5
  which_sheet <- ifelse(plate <= 4, 'plates1_4', 'plate5')
  # Import the sheet
  plate_info[[plate]] <- read_sheet(
    sheets_list[[which_sheet]],
    sheet = paste0("plate", plate),
    range = paste0('plate', plate, '!A3:P99')
    ) %>%
    mutate(
      'plate' = paste0('plate', plate)
    )
}
plate_info <- do.call('rbind', plate_info)


plate_info <- plate_info %>%
  rename_with(~ gsub(' ', '_', .x)) %>%
  mutate(
    Sample_name = gsub(' ', '_', Sample_name),
    barcode = paste0(Barcode7, Barcode5),
    NGS_sample_ID = as.character(NGS_sample_ID)
    ) %>%
  left_join(., bam_info, by = c('NGS_sample_ID', 'barcode'))

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
  to   = paste0(plate_info$directory, plate_info$Sample_name, ".bam")
)
file.rename(
  from = paste0(plate_info$directory, plate_info$filename, '.bai'),
  to   = paste0(plate_info$directory, plate_info$Sample_name, ".bam.bai")
)


write_csv(plate_info, "03_processing/plate_info.csv")
