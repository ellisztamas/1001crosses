#' Exploratory plots for the parental SNP matrix
#'
#' This plots the results of 03_processing/01_parental_SNP_matrix/01_inspect_initial_snp_matrix.sh
#'
#' Tom Ellis, adpating code from https://speciationgenomics.github.io/filtering_vcfs/
#' 3rd Janurary 2024

# load tidyverse package
library(tidyverse)

# Directory with summary tables.
data_dir <- "03_processing/01_parental_SNP_matrix/summary_stats_initial"

# Variant quality
# Ranges from 30 to 3e6
var_qual <- read_delim(
  paste0(data_dir, "/data_inspection.lqual"),
  delim = "\t",
  col_names = c("chr", "pos", "qual"), skip = 1)
var_qual %>%
  ggplot(aes(qual)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()

# Variant depth ranges from ~0 to ~68, with a peak at 19
var_depth <- read_delim(
  paste0(data_dir, "/data_inspection.ldepth.mean"),
  delim = "\t",
  col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
var_depth %>%
  ggplot(aes(mean_depth)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()

# Allele frequencies
# These are strongly skewed left, with a mean of 3%. That is weird to be honest.
# They are likley to be singletons.
var_freq <- read_delim(
  paste0(data_dir, "/data_inspection.frq"),
  delim = "\t",
  col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq$maf <- var_freq %>%
  select(a1, a2) %>%
  apply(1, function(z) min(z))
var_freq %>%
  # filter(maf < 0.1) %>%
  ggplot(aes(maf)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
summary(var_freq$maf)

# Most loci have only a few percent missing data, with a long tail.
# It should be fine to filter at 10%
var_miss <- read_delim(
  paste0(data_dir, "/data_inspection.lmiss"),
  delim = "\t",
  col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
var_miss %>%
  ggplot(aes(fmiss)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
summary(var_miss$fmiss)

# Heterozygosity per variant
# Mostly close to zero, with a long tail up to 1.
var_het <- read_delim(
  paste0(data_dir, "/heterozygosity_per_site.tsv"),
  delim = "\t")
var_het %>%
  ggplot(aes(het_rate)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
# Create a table of sites with heterozygosity >5%
var_het %>%
  filter(het_rate > 0.01) %>%
  select(CHROM, POS) %>%
  write_tsv(
    "03_processing/01_parental_SNP_matrix/data_inspection/heterozygous_sites_to_purge.tsv",
    col_names = FALSE)

# Depth across individuals
ind_depth <- read_delim(
  paste0(data_dir, "/data_inspection.idepth"),
  delim = "\t",
  col_names = c("ind", "nsites", "depth"), skip = 1)
ind_depth %>%
  ggplot(aes(depth)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()

# Missing data across individuals
# Peak at ~6%, with a tail up to 0.37
ind_miss  <- read_delim(
  paste0(data_dir, "/data_inspection.imiss"),
  delim = "\t",
  col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1
  )
ind_miss %>%
  ggplot(aes(fmiss)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()

# Heterozygosity across individuals
# Most are >0.9, with a tail down to 0.37
ind_het <- read_delim(
  paste0(data_dir, "/data_inspection.het"),
  delim = "\t",
  col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
ind_het %>%
  ggplot(aes(f)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()


