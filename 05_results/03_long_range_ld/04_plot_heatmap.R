#' Plot long-range LD in the parents and progeny
#'
#' Creates a matrix of positions within and between chromosomes as
#' a facet grid and plots positions with substantial LD in the
#' parents and all F8 genotypes.
#'
#' Tom Ellis, 5th May 2025

library(tidyverse)

# CSV files for the parents and F8s giving Chr, position and r2
# This is for loci with a minimum r2; see previous scripts in this directory
parents <- read_csv(
  "05_results/03_long_range_ld/output/parents_snps_in_LD.csv",
  col_types = 'ccd'
  ) %>%
  mutate(
    generation = "parents"
  )
progeny <- read_csv(
  "05_results/03_long_range_ld/output/progeny_snps_in_LD.csv",
  col_types = 'ccd'
  ) %>%
  mutate(
    generation = "progeny"
  )

# Join LD datasets so they have a common coordindate system.
r2_table <- rbind(parents, progeny) %>%
  separate(i, into=c("chr_i", 'pos_i'), sep=":") %>%
  separate(j, into=c("chr_j", 'pos_j'), sep=":")

# Set a minimum value of r2 to include in the plot.
minimum_r2_to_plot <- 0.5

r2_table %>%
  # Set a minimum r2 to plot
  filter(r2 > minimum_r2_to_plot | r2 > minimum_r2_to_plot) %>%
  # To get the values for the progeny to plot in the lower-right
  # corner, swap the x and y coordinates (pos_i and pos_j).
  # This is done by creating temporary variables chr_i_tmp and pos_i_tmp
  mutate(
    chr_i_tmp = chr_i,
    pos_i_tmp = pos_i
  ) %>%
  # Now swap the values of pos_i and pos_j, and chr_i and chr_j
  # for the progeny generation
  mutate(
    pos_i = ifelse(generation == "progeny", pos_j,     pos_i),
    pos_j = ifelse(generation == "progeny", pos_i_tmp, pos_j),
    chr_i = ifelse(generation == "progeny", chr_j,     chr_i),
    chr_j = ifelse(generation == "progeny", chr_i_tmp, chr_j),
  ) %>%
  # Make base positions numeric and on a megabase scale
  mutate(
    pos_i = as.numeric(pos_i)/1e6,
    pos_j = as.numeric(pos_j)/1e6,
  ) %>%
  # Plot a 5x5 grid of LD across the genome
  ggplot(aes(pos_i, pos_j, colour=r2)) +
  geom_point(size=0.3) +
  labs(
    x = "Position (Mb)",
    y = "Position (Mb)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  scale_color_binned(type="viridis", direction=-1, n.breaks = 6) +
  facet_grid(fct_rev(chr_j) ~ chr_i, scales="free")

ggsave(
  filename = "05_results/03_long_range_ld/output/long_range_LD.png",
  device = "png",
  height = 16.9, width = 16.9, units = "cm"
  )


# Histograms of LD within and between chromosomes
r2_table %>%
  mutate(
    within_between = ifelse(chr_i == chr_j, "within", "between")
  ) %>%
  # filter(r2 > 0.1) %>%
  ggplot(aes(x=r2, colour=generation)) +
  geom_freqpoly() +
  theme_bw() +
  facet_grid(~within_between)

# Scatter plot of between-chromosome LD
r2_table %>%
  pivot_wider(names_from = generation, values_from = r2) %>%
  # filter(chr_i != chr_j) %>%
  mutate(
    within_between = ifelse(chr_i == chr_j, "within", "between")
  ) %>%
  ggplot(aes(x=parents, y= progeny, colour = within_between)) +
  geom_point()
