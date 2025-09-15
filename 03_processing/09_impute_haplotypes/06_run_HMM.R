#!/usr/bin/env Rscript

library("vcfR")
library("optparse")
library("readr")
library("ggplot2")
# Functions to run the HMM
source("03_processing/09_impute_haplotypes/04_emission_probabilities.R")
source("03_processing/09_impute_haplotypes/05_viterbi_hmm.R")

# Import command-line arguments
option_list <- list(
  make_option(c("-v", "--vcf_file"),
              type="character",
              help="Path to a VCF file, containing both parents and the progeny."),
  make_option(c("-m", "--parent1"),
              type="character",
              help="Name of the first parent."),
  make_option(c("-f", "--parent2"),
              type="character",
              help="Name of the second parent."),
  make_option(c("-x", "--progeny"),
              type="character",
              help="Name of the progeny."),
  make_option(c("-e", "--mapping_error"),
              type="numeric",
              help="Probability that a read aligns to the wrong place."),
  make_option(c("-o", "--outdir"),
              type="character",
              help="Output directory to write the output files.")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Read VCF file
vcf <- read.vcfR(opt$vcf_file)

# Calculate emission probabilities at each marker that the progeny individual is
# homozygous for either parental acnestry, or heterozygous.
emission_data <- emission_probabilities(
  vcf_file = vcf,
  parent1 = opt$parent1,
  parent2 = opt$parent2,
  progeny = opt$progeny,
  err = opt$mapping_error
  )
# Run the Viterbi algoirthm on the data.
# `remove_singletons` purges haplotypes that only include a single marker,
# because I cannot conveive that these are real.
hidden_states <- viterbi_hmm(emission_data, remove_singletons = TRUE)

# Write tables to disk
write_tsv(
  hidden_states$markers,
  file = paste0(opt$outdir, "/", opt$progeny, "_marker_states.tsv")
)
write_tsv(
  hidden_states$haplotypes,
  file = paste0(opt$outdir, "/", opt$progeny, "_haplotypes.bed")
)


