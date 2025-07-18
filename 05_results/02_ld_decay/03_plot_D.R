# Summarise LD within bins.

# Tom Ellis, adapting from code by Jinliang Yang:
# https://jyanglab.com/AGRO-932/chapters/a2.1-qg/rex11_gwas2.html#12

library("data.table")
library(plyr)
library(tidyverse)
# library(ggplot2)

ld_files <- list(
  parents="05_results/02_ld_decay/output/parental_lines_D.ld.gz",
  progeny='05_results/02_ld_decay/output/F8_phased_imputed_D.ld.gz'
)



# Window size to average over
bin_size = 100
ld_bins <- vector('list', 3)
names(ld_bins) <- names(ld_files)
for(name in names(ld_files)){
  # Get mean r2 within 100bp windows.
  filename <- ld_files[[name]]
  df <- fread(filename, data.table = FALSE) %>%
    mutate(
      dist = BP_B - BP_A,
      bin  = round(dist / bin_size, 0)
    )

  ld_bins[[name]] <- ddply(
    df, .(bin), summarise,
    meanD = mean(D)
  ) %>%
    mutate(
      Generation = name,
      distance = bin*100 +50
    )
}



do.call(what = 'rbind', ld_bins) %>%
  ggplot(aes(x=distance/1000, y = meanD, colour=Generation)) +
  geom_point(size=0.5) +
  labs(
    x = "Physical distance (kb)",
    y = "Mean D"
  ) +
  theme_bw()


do.call(what = 'rbind', ld_bins) %>%
  pivot_wider(names_from = Generation, values_from = meanD) %>%
  mutate(
    rfrac_3 = (distance*1e-6 * 5) / 10,
    expD_3  = (1-rfrac_3)*(1-(2*rfrac_3))^7 * parents
  ) %>%
  select(-rfrac_3) %>%
  pivot_longer(parents:progeny, names_to = 'Generation', values_to = "meanD") %>%
  ggplot(aes(x=distance/1000, y = meanD, colour=Generation)) +
  geom_point(size=0.5) +
  geom_line(aes(x = distance/1000, y=expD_3), colour='black') +
  labs(
    x = "Physical distance (kb)",
    y = "Mean D"
  ) +
  theme_bw()


# ggsave(
#   "05_results/02_ld_decay/output/short_range_LD.png",
#   device = "png",
#   units = "in",
#   height = 6, width =8
# )
