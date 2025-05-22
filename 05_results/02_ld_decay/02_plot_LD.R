# Summarise LD within bins.

# Tom Ellis, adapting from code by Jinliang Yang:
# https://jyanglab.com/AGRO-932/chapters/a2.1-qg/rex11_gwas2.html#12

library("data.table")
library(plyr)
library(tidyverse)
# library(ggplot2)

ld_files <- list(
  parents='05_results/02_ld_decay/output/parental_lines.ld.gz',
  progeny="05_results/02_ld_decay/output/F8_phased_imputed.ld.gz"
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
    meanr2 = mean(R2)
  ) %>%
    mutate(
      Generation = name
    )
}

# ld_bins <- lapply(names(ld_files), function(name){
#   filename <- ld_files[[name]]
#   fread(filename, data.table = FALSE) %>%
#     mutate(
#       dist = BP_B - BP_A,
#       bin  = round(dist / bin_size, 0)
#     ) %>%
#     plyr::ddply(
#       .(bin), summarise, mean_r2 = mean(R2)
#     ) %>%
#     mutate(
#       dataset = name
#     )
# })

do.call(what = 'rbind', ld_bins) %>%
  ggplot(aes(x=bin/10, y = meanr2, colour=Generation)) +
  geom_point() +
  labs(
    x = "Physical distance (kb)",
    y = expression(paste('Mean ', r^{2}))
  ) +
  lims(
    # x=c(0,25),
    y=c(0,0.6)
    ) +
  theme_bw()

ggsave(
  "05_results/02_ld_decay/output/short_range_LD.png",
  device = "png",
  units = "in",
  height = 6, width =8
  )
