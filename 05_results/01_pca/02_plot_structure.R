library("tidyverse")

library('ggpubr')

# List of file paths giving eigenvectors
eigenvec_files <- list(
  rep1 = "05_results/01_pca/output/F8_snp_matrix_purged_rep1.eigenvec",
  rep2 = "05_results/01_pca/output/F8_snp_matrix_purged_rep2.eigenvec",
  parents = "05_results/01_pca/output/parental_snp_matrix.eigenvec"
)
# List of file paths giving eigen values
eigenval_files <- list(
  rep1 = "05_results/01_pca/output/F8_snp_matrix_purged_rep1.eigenval",
  rep2 = "05_results/01_pca/output/F8_snp_matrix_purged_rep2.eigenval",
  parents = "05_results/01_pca/output/parental_snp_matrix.eigenval"
)

# Import PCA data
eigenvecs <- lapply(names(eigenvec_files), function(name){
  filename <- eigenvec_files[[name]]
  read.table(filename, header=TRUE) %>%
    as_tibble() %>%
    mutate(
      dataset = name
    )
})
eigenvals <- lapply(names(eigenval_files), function(name){
  filename <- eigenval_files[[name]]
  read.table(filename, col.names="eigenval", header=FALSE) %>%
    as_tibble() %>%
    mutate(
      eigenval = round(100*(eigenval / sum(eigenval)), 1),
      dataset = name
    )
})

ggarrange(
  eigenvecs[[3]] %>%
    ggplot(aes(x=PC1, y = PC2)) +
    geom_point() +
    labs(
      x = paste0("PC1 (", eigenvals[[3]]$eigenval[1], "%)"),
      y = paste0("PC2 (", eigenvals[[3]]$eigenval[2], "%)"),
      title = "Parents"
    ) +
    theme_bw(),

  eigenvecs[[1]] %>%
    ggplot(aes(x=-PC1, y = PC2)) +
    geom_point() +
    labs(
      x = paste0("PC1 (", eigenvals[[1]]$eigenval[1], "%)"),
      y = paste0("PC2 (", eigenvals[[1]]$eigenval[2], "%)"),
      title = "F8 (rep. 1)"
    ) +
    theme_bw(),

  eigenvecs[[2]] %>%
    ggplot(aes(x=-PC1, y = -PC2)) +
    geom_point() +
    labs(
      x = paste0("PC1 (", eigenvals[[2]]$eigenval[1], "%)"),
      y = paste0("PC2 (", eigenvals[[2]]$eigenval[2], "%)"),
      title = "F8 (rep. 2)"
    ) +
    theme_bw(),

  ncol =3
)

ggsave(
  filename = "05_results/01_pca/output/population_structure.png",
  device = "png",
  units = 'cm', height = 8, width = 16.9
)

