#' Script to plot the distribution of LD stats between chromosomes.
#' It doesn't make sense to use these data for within-chomosome comparisons
#' because the loci are filtered to be >10kb apart
#' 
#' Tom Ellis 18th July 2025

library('tidyverse')

# Load text files giving LD between pairs of SNPs genome wide.
# They are massive, so the `sed` command first filters these to leave only every
# 10000th pairs
parents <- read_csv(
  pipe(
    "sed -n '1~10000p' 05_results/03_long_range_ld/output/parents_snps_in_LD.csv"
  ),
  col_names = c("i", "j", "d", "dprime", "r2"),
  col_types = 'ccddd'
) %>% 
  mutate(generation = "parents")

progeny <- read_csv(
  pipe(
    "sed -n '2~10000p' 05_results/03_long_range_ld/output/progeny_snps_in_LD.csv"
  ),
  col_names = c("i", "j", "d", "dprime", "r2"),
  col_types = 'ccddd'
) %>% 
  mutate(generation = "progeny")

# Join LD datasets so they have a common coordindate system.
r2_table <- rbind(parents, progeny) %>%
  separate(i, into=c("chr_i", 'pos_i'), sep=":") %>%
  separate(j, into=c("chr_j", 'pos_j'), sep=":")

# Histograms of r2 between chromosomes
r2_table %>%
  filter(chr_i != chr_j) %>% 
  ggplot(aes(x=r2, colour=generation)) +
  geom_freqpoly() +
  lims(
    x=c(0,1)
  ) +
  labs(
    title = "Between chromosomes",
    x = expression(paste('Linkage disequliibrium (', r^{2}, ")")),
    y = "Pairs of loci"
  ) +
  theme_bw()

# D between chromosomes
r2_table %>%
  filter(chr_i != chr_j) %>% 
  ggplot(aes(x=d, colour=generation)) +
  geom_freqpoly() +
  labs(
    title = "Between chromosomes",
    x = "D",
    y = "Pairs of loci"
  ) +
  theme_bw()
