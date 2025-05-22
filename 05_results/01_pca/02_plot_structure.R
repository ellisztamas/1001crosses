library("tidyverse")

library('ggpubr')

# List of file paths giving eigenvectors
eigenvec_files <- list(
  parents = "05_results/01_pca/output/parental_lines.eigenvec",
  progeny = "05_results/01_pca/output/F8_phased_imputed.eigenvec"
)
# List of file paths giving eigen values
eigenval_files <- list(
  parents = "05_results/01_pca/output/parental_lines.eigenval",
  progeny = "05_results/01_pca/output/F8_phased_imputed.eigenval"
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
names(eigenvecs) <- names(eigenvec_files)

eigenvals <- lapply(names(eigenval_files), function(name){
  filename <- eigenval_files[[name]]
  read.table(filename, col.names="eigenval", header=FALSE) %>%
    as_tibble() %>%
    mutate(
      eigenval = round(100*(eigenval / sum(eigenval)), 1),
      dataset = name
    )
})
names(eigenvals) <- names(eigenval_files)

ggarrange(
  eigenvecs$parents %>%
    ggplot(aes(x=PC1, y = PC2)) +
    geom_point() +
    labs(
      x = paste0("PC1 (", eigenvals$parents$eigenval[1], "%)"),
      y = paste0("PC2 (", eigenvals$parents$eigenval[2], "%)"),
      title = "Parents"
    ) +
    theme_bw(),

  eigenvecs$progeny %>%
    ggplot(aes(x=PC1, y = PC2)) +
    geom_point() +
    labs(
      x = paste0("PC1 (", eigenvals$progeny$eigenval[1], "%)"),
      y = paste0("PC2 (", eigenvals$progeny$eigenval[2], "%)"),
      title = "Progeny"
    ) +
    theme_bw(),

  ncol =2
)

ggsave(
  filename = "05_results/01_pca/output/population_structure.png",
  device = "png",
  units = 'cm', height = 8, width = 16.9
)



