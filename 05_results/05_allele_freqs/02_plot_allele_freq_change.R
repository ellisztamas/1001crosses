


library("data.table")
library(tidyverse)
library(ggpubr)

frq_paths <- list(
  rep1    = "05_results/05_allele_freqs/output/F8_snp_matrix_purged_rep1.afreq",
  rep2    = "05_results/05_allele_freqs/output/F8_snp_matrix_purged_rep2.afreq",
  parents = "05_results/05_allele_freqs/output/parental_snp_matrix.afreq"
)

frq_list <- lapply(frq_paths, function(path){
  fread(path, data.table = FALSE) %>%
    dplyr::select(`#CHROM`, ID, ALT_FREQS) %>%
    rename(
      chr = `#CHROM`,
      id = ID,
      freq = ALT_FREQS
    )
}
)

allele_freqs <- frq_list$rep1 %>%
  inner_join(
    frq_list$rep2, by = c("chr", 'id'), suffix = c("_rep1", '_rep2')
  ) %>%
  inner_join(
    frq_list$parents, by = c("chr", 'id')
  ) %>%
  rename(
    freq_parents = freq
  )

ggarrange(
  allele_freqs %>%
    ggplot(aes(x = freq_parents, y = freq_rep1, colour = chr)) +
    geom_point(),

  allele_freqs %>%
    ggplot(aes(x = freq_parents, y = freq_rep2, colour = chr)) +
    geom_point()
)

ggsave(
  filename = "05_results/05_allele_freqs/output/allele_freqs.png",
  device = "png",
  units = 'cm', height = 8, width = 16.9
)
