#' Select a set of SNPs to simulate GWAS results
#'
#' Input:
#'    Table summarising allele frequency and correlations with PCs of population
#'      structure for each SNP. Created by 01_allele_correlations.py.
#' Output:
#'    TSV file giving SNPs for further simulation. This is a subset of the input
#'      file.
#'
#' The choice of SNPs might change, but at present this bins SNPs by allele freq
#' and the maximum correlation with any PC of population structure, and chooses
#' up to some number from each bin. I do not currently check whether things are
#' linked, but the input matrix is sufficiently sparse that I don't need to.
#'
#' Tom Ellis, 10th July 2024

library(tidyverse)

# Load the data
variant_table <- read_csv(
  "05_results/09_gwas_single_loci/output/variant_table.csv", show_col_types = FALSE
  )

# Tody up
variant_table <- variant_table %>%
  # Filter for missing data
  filter(
    complete.cases(.),
    maf > 0,
    mac >= 400
  ) %>%
  # Divide the SNPs into bins allele-frequency and correlation with any PC.
  # Notice that the first bin for correlations is much bigger than the others
  mutate(
    maf = ifelse(maf > 0.5, 1-maf, maf),
    af_bin = cut(maf, c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)),
    max_cor_bin = cut(max_cor, c(0,0.5, 0.6, 0.7, 0.8, 0.9, 1)),
    mean_cor_bin = cut(mean_cor, seq(0,0.5, 0.1))
  ) %>%
  filter(
    !is.na(af_bin), !is.na(max_cor_bin)
  )

# Plot the number of SNPs in each max_cor bin
variant_table %>%
  group_by(af_bin, max_cor_bin) %>%
  summarise(
    n = n()
  ) %>%
  ggplot(aes(x=af_bin, y = n)) +
  geom_point() +
  facet_grid(~max_cor_bin)
# Now by mean correlation with each PC
variant_table %>%
  group_by(af_bin, mean_cor_bin) %>%
  summarise(
    n = n()
  ) %>%
  ggplot(aes(x=af_bin, y = n)) +
  geom_point() +
  facet_grid(~mean_cor_bin)
# Mean and maximum correlation are themselves correlated, and this doesn't
# depend on MAF, so stick with max
variant_table %>%
  ggplot(aes(x=mean_cor, y = max_cor)) +
  geom_point() +
  facet_grid(!af_bin)


# Select SNPs from each bin.
# Randomly choose up to this many:
max_sample_size <- 20
# If there aren't that many in a bin, use all that there are.

# Choose SNPs
set.seed(215)
snps_to_simulate <- variant_table %>%
  group_split(max_cor_bin, af_bin) %>%
  lapply(., function(x){
    n <- ifelse(nrow(x) >= max_sample_size, max_sample_size, nrow(x))
    x[sample(1:nrow(x), size=n, replace = FALSE),]
    }) %>%
  bind_rows() %>%
  arrange(chrom, pos)

# Write to disk
snps_to_simulate %>%
  write_delim(
    "05_results/09_gwas_single_loci/output/snps_to_simulate.tsv",
    delim = "\t"
    )
