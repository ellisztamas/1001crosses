# Summarise LD within bins.

# Tom Ellis, adapting from code by Jinliang Yang:
# https://jyanglab.com/AGRO-932/chapters/a2.1-qg/rex11_gwas2.html#12

library("data.table")
library(tidyverse)

ld_files <- list(
  rep1="05_results/02_ld_decay/output/F8_snp_matrix_purged_rep1.ld.gz",
  rep2="05_results/02_ld_decay/output/F8_snp_matrix_purged_rep2.ld.gz",
  parents='05_results/02_ld_decay/output/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.ld.gz'
)

# Window size to average over
bin_size = 100
# Get mean r2 within 100bp windows.
ld_bins <- lapply(names(ld_files), function(name){
  filename <- ld_files[[name]]
  fread(filename, data.table = FALSE) %>%
    mutate(
      dist = BP_B - BP_A,
      bin  = round(dist / bin_size, 0)
    ) %>%
    ddply(
      .(bin), summarise, mean_r2 = mean(R2)
    ) %>%
    mutate(
      dataset = name
    )
})

do.call(what = 'rbind', ld_bins) %>%
  ggplot(aes(x=bin*100, y = mean_r2, colour=dataset)) +
  geom_point() +
  labs(
    x = "Physical distance (bp)",
    y = "Mean R2"
  ) +
  lims(x=c(0,25000))

