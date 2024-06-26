#' Line up original and corrected sample names of the F8s.
#'
#' Creates a dataframe giving original sample names of the F8s with those
#' after Pieter checked them. This also filters out those Pieter deemed
#' dubious.
#'
#' Sample names are merged based on the BAM filenames they are associated with.
#'
#' Inputs:
#'    Pieter's original sample sheet
#'    Pieter's sample sheet after he had checked validated geneotypes, including
#'      updated sample names and a column "genoSelection" stating whether the
#'      sample could be validated.
#' Outputs:
#'    A data frame with a column for original names and a column for updated
#'    names. This is not saved to disk.
#'
#' Tom Ellis 25th June 2024

library(tidyverse)

original_plate_info <- read_csv(
  "../crosses_continued/004.F8/001.genotyping/001.data/sample.info.csv"
  )%>%
  select(Sample_name, filename) %>%
  filter(!is.na(filename))


corrected_plate_info <- read_csv(
  "../crosses_continued/004.F8/002.mapping/001.data/plate_info.csv"
) %>%
  filter(
    !is.na(filename),
    genoSelection == "yes"
    ) %>%
  select(Sample_name, filename)

sample_name_corrections <- inner_join(
  original_plate_info, corrected_plate_info,
  by = "filename",
  suffix = c("_original", '_corrected')
  ) %>%
  mutate(
    Sample_name_original =  str_replace(Sample_name_original, "_F8_", "_"),
    Sample_name_corrected = str_replace(Sample_name_corrected, "_F8_", "_")
  ) %>%
  select(
    -filename
  )

