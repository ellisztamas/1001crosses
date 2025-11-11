#!/usr/bin/env Rscript

#' Script to calculate the inflation statistics on a GWAS results file.
#'
#' This takes a GEMMA results files (ending in `.assoc.txt`) and calculates
#' lambda on the column 'p_lrt'. It also counts the number of 500kb windows
#' that contain at least one Bonferroni-corrected significant association.
#'
#' The output is a dataframe with a single row giving cohort, replicate ID,
#' number of loci, (known) heritability, and the two inflation statistics.
#'
#' Tom Ellis 4th November 2025


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments provided
if (length(args) != 4) {
  cat("Usage: Rscript script.R --in input_file --out output_file\n")
  cat("   or: Rscript script.R input_file output_file\n")
  quit(status = 1)
}

# Parse arguments
if (args[1] == "--in" && args[3] == "--out") {
  # Named arguments: --in file1 --out file2
  input_file <- args[2]
  output_file <- args[4]
} else if (length(args) == 2) {
  # Positional arguments: file1 file2
  input_file <- args[1]
  output_file <- args[2]
} else {
  cat("Error: Invalid arguments\n")
  quit(status = 1)
}



# Main  -------------------------------------------------------------------

suppressPackageStartupMessages(
  library(tidyverse)
)
library(filelock)

gwas <- read_tsv(input_file, show_col_types = FALSE)

# Prepare a tibble with cohort, replicate ID, number of loci and heritability
params <- basename(input_file) %>%
  gsub(".assoc.txt", "", .) %>%
  str_split_fixed("_", 4)

out <- tibble(
  cohort    = params[1],
  rep       = params[2],
  maf_range = params[3],
  h2        = params[4]
)

# Calculate lambda
chisq_stats <- qchisq(1 - gwas$p_lrt, df = 1)
median_chisq <- median(chisq_stats, na.rm = TRUE)
expected_median <- qchisq(0.5, df = 1)
lambda_gc <- median_chisq / expected_median

out$lambda <- round(lambda_gc, 4)


# Count how many windows have a significant result
# Create a column giving 500kb windows for each chromosome.
gwas$window <- cut(gwas$ps, seq(1, max(gwas$ps), 500000))
gwas$window <- paste0(gwas$chr, gwas$window)
# How many windows have at least one SNP with p_lrt < Bonferoni threshold.
bonf_threshold <- 0.05/nrow(gwas)
n_significant_windows <- length(
  unique(
    gwas$window[gwas$p_lrt < bonf_threshold]
    )
  )

out$n_windows <- n_significant_windows


# Write the output
# This needs to use a filelock to stop lines getting mangled.
lock_path <- paste0(output_file, ".lock")
# Acquire the lock (waits until available)
lock <- lock(lock_path)

# Now safe to write
out %>%
  write_tsv(
    file = output_file,
    append=TRUE
    )
# Release the lock
unlock(lock)
