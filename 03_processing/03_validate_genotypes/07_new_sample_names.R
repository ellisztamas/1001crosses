#!/usr/bin/env Rscript
#
#' Create the text files needed to rename VCF files.
#'
#' The VCF file created in 03_processing/02_original_sample_sheet gives sample
#' names as their sequencing plate, and position on those plates (eg 2_G3).
#' This script takes a text file containing those names, and returns:
#'     1. A text file of the same length with sample names (eg 6201x9399 F8 rep1)
#'        in the same order
#'     2. Separate files for sameples from reps 1 and 2, excluding samples
#'        found to have no data, be self-fertilised, or be the wrong species.

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please supply an argument giving the working directory", call.=FALSE)
}
workdir <- args[1]
# workdir <- "03_processing/03_validate_genotypes/output"

suppressMessages(
  library('tidyverse', verbose = FALSE)
)

# Text file giving sample names as they stand in the VCF file as positions
# e.g. 4_G3
position_names <- read_tsv(
  paste0(workdir, "/vcf_sample_names_as_positions.txt"),
  col_names = "pos_name"
  )
# Results of validating the genotypes
ibdresults <- read_csv(
  "03_processing/03_validate_genotypes/output/ibdpainting_results.csv"
  ) %>%
  mutate(
    pos_name= paste0(plate, "_", row, col),
    new_name = gsub(" F8 ", "_", new_name)
  )

# Merge
corrected_names <- position_names %>%
  left_join(ibdresults, by='pos_name')

# Text file giving sample names in the order they appear in the VCF file
x <- corrected_names %>%
  select(new_name) #%>%
  write_csv(paste0(workdir, "/vcf_sample_names_as_names.txt"), col_names = FALSE)

# Remove samples that had no data, were self-fertilised, or were the wrong species
y <- corrected_names %>%
  filter(
    diagnosis != "Self-pollinated",
    ! grepl("Aa1xAa2", name),
    diagnosis != "no data"
  ) %>%
  select(new_name) #%>%
  write_csv(paste0(workdir, "/vcf_sample_names_filtered.txt"), col_names = FALSE)

y$new_name %in%  x$new_name
