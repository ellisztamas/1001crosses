#' Principle component decomposition of the QTL region on Chr5
#'
#' Import the matrix of SNP genotypes in the region Chr5:16207115-16248614.
#' Run a PCA on that, and save PC1 as a text file that GEMMA can read as a
#' table of covariates.
#'
#' Pay attention to the hacky imputation of missing data.
#'
#' Tom Ellis 16th October 2025

qtl_region <- read_csv(
  "05_results/23_admixed_only_GWAS/output/05_extract_QTL_region/Chr5_QTL_region.txt",
  show_col_types = FALSE
)

# Matrix of genotypes
# Rows as lines, columns as SNPs
mat <- qtl_region %>%
  pivot_longer(`1002x6244_rep1`:`992x9332_rep2`) %>%
  pivot_wider(names_from = snp) %>%
  tibble::column_to_rownames(var = 'name') %>%
  as.matrix()

# There are NAs in the matrix.
# I think most of them should be 0s.
# A very hacky solution is to SNPs with more than half missing data,
# and set remaining NAs to the SNP-mean value.
# Please don't tell Richard McElreath I did this.
colMeans(is.na(mat)) > 0.5
mat_imputed <- apply(mat[,colMeans(is.na(mat)) < 0.3], 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})

# PCA on the imputed matrix.
pca_qtl_region <- prcomp(mat_imputed)

# Tibble that can be used as a covariate file by GEMMA.
# Column of 1s, line names, and PC1
tibble(
  line = row.names(pca_qtl_region$x),
  PC1 = pca_qtl_region$x[,1]
) %>%
  right_join(pheno, by ='line') %>%
  select(FID, PC1) %>%
  mutate(
    FID=1
  ) %>%
  write_tsv(
    "05_results/23_admixed_only_GWAS/output/Chr5_QTL_region_PC1.tsv",
    col_names = FALSE
    )


