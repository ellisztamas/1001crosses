#' An attempt to retrieve allele frequencies from GWAS results


gwas.cross <- read_delim('05_results/07_gemma_seed_size/output/with_K/seed_size_blups_F9_rep1.assoc.txt', delim = '\t')
gwas.parent <- read_delim('05_results/07_gemma_seed_size/output/with_K/seed_size_blups_parents.assoc.txt', delim = '\t')

gwas.result <- inner_join(gwas.cross, gwas.parent, suffix = c('.cross', '.parent'), by = c('chr', 'ps'))

maf.plt <- gwas.result %>%
  ggplot(aes(x = af.parent, y = af.cross)) +
  geom_point(size = 0.1) +
  theme_classic() +
  labs(
    x = 'parents',
    y = 'crosses',
    title = 'minor allele frequencies'
    ) +
  geom_abline(intercept = 0, slope = 1, col = 'red')

ggsave(plot = maf.plt,
       filename = '05_results/05_allele_freqs/output/af_from_gwas.png',
       type = 'cairo')
