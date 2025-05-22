library(tidyverse)
library(ggpubr)

# Vector of GEMMA output files.
gemma_output_paths <- Sys.glob(
  "05_results/16_gemma_ft12/output/no_K/*.assoc.txt"
)
# Drop the file with all F9s combined.
gemma_output_paths <- grep(
  "F9_combined", gemma_output_paths,
  invert = TRUE, value = TRUE
)


gwas <- vector("list", length=length(gemma_output_paths))
names(gwas) <- c("Parents", "F9_cohort_1", "F9_cohort_2")

for(f in 1:length(gwas)){
  gemma_file <- read_tsv(
    gemma_output_paths[f],
    # n_max = 1000,
    col_types = 'cciiiccnnnnnn'
  ) %>%
    mutate(
      log10p = -log10(p_lrt),
      varexp = beta^2 * af * (1-af)
    )

  gwas[[f]] <- gemma_file

}

# Merge the result files
merge_gwas <- inner_join(
  gwas$Parents, gwas$F9_cohort_1,
  by = c("chr", "ps"),
  suffix = c("_parents", "")
) %>%
  inner_join(
    gwas$F9_cohort_2,
    by = c("chr", "ps"),
    suffix = c("_cohort1", "_cohort2")
  )


# Spearman rank correlations for p-value distributions
pval_correlations <- c(
  cor(merge_gwas$log10p_parents, merge_gwas$log10p_cohort1, method='s'),
  cor(merge_gwas$log10p_parents, merge_gwas$log10p_cohort2, method='s'),
  cor(merge_gwas$log10p_cohort1, merge_gwas$log10p_cohort2, method='s')
)
pval_correlations <- sapply(pval_correlations, round, 3)


ggpubr::ggarrange(
  merge_gwas %>%
    ggplot(aes(x = log10p_parents, y=log10p_cohort1 )) +
    geom_point()+
    labs(
      x = "-log10 p (parents)",
      y = "-log10 p (F9 cohort 1)",
      title = "Parents vs F9 cohort 1",
      subtitle = paste("\u03c1 = ", pval_correlations[1])
    ),
  merge_gwas %>%
    ggplot(aes(x = log10p_parents, y=log10p_cohort2 )) +
    geom_point()+
    labs(
      x = "-log10 p (parents)",
      y = "-log10 p (F9 cohort 2)",
      title = "Parents vs F9 cohort 2",
      subtitle = paste("\u03c1 = ", pval_correlations[2])
    ),
  merge_gwas %>%
    ggplot(aes(x = log10p_cohort1, y=log10p_cohort2 )) +
    geom_point()+
    labs(
      x = "-log10 p (F9 cohort 1)",
      y = "-log10 p (F9 cohort 2)",
      title = "Cohorts 1 and 2",
      subtitle = paste("\u03c1 = ", pval_correlations[3])
    ),

  nrow=1, ncol=3
)

ggsave(
  filename = "05_results/16_gemma_ft12/output/compare_pvals.png",
  device = "png",
  units = "cm",
  height  = 15, width = 16.9
)


# Spearman rank correlations for effect-size distributions
beta_correlations <- c(
  cor(merge_gwas$beta_parents, merge_gwas$beta_cohort1, method='s'),
  cor(merge_gwas$beta_parents, merge_gwas$beta_cohort2, method='s'),
  cor(merge_gwas$beta_cohort1, merge_gwas$beta_cohort2, method='s')
)
beta_correlations <- sapply(beta_correlations, round, 3)



ggpubr::ggarrange(
  merge_gwas %>%
    ggplot(aes(x = beta_parents, y=beta_cohort1 )) +
    geom_point()+
    labs(
      x = "SNP effect (parents)",
      y = "SNP effect (F9 cohort 1)",
      title = "Parents vs F9 cohort 1",
      subtitle = paste("\u03c1 = ", beta_correlations[1])
    ),

  merge_gwas %>%
    ggplot(aes(x = beta_parents, y=beta_cohort2 )) +
    geom_point()+
    labs(
      x = "SNP effect (parents)",
      y = "SNP effect (F9 cohort 2)",
      title = "Parents vs F9 cohort 2",
      subtitle = paste("\u03c1 = ", beta_correlations[2])
    ),

  merge_gwas %>%
    ggplot(aes(x = beta_cohort1, y=beta_cohort2 )) +
    geom_point()+
    labs(
      x = "SNP effect (F9 cohort 1)",
      y = "SNP effect (F9 cohort 2)",
      title = "Cohorts 1 and 2",
      subtitle = paste("\u03c1 = ", beta_correlations[3])
    ),

  nrow=1, ncol=3
)

ggsave(
  filename = "05_results/16_gemma_ft12/output/compare_betas.png",
  device = "png",
  units = "cm", height  = 15, width = 16.9
)



