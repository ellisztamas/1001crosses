#' Simulate phenotypes based on observed genotypes.
#'
#' This script imports a matrices of SNP genotype calls for the parents
#' and two cohorts of F8s. It is assumed the SNPs are the same for the three
#' cohorts.
#' To generate phenotypes the script:
#'   - Transposes that matrix
#'   - Centres and standardises variants at each marker.
#'   - Gives each marker is an effect size, which is multiplied by standardised
#'       genotype values, and summed across markers (in this case there is only
#'       one marker per phenotype, but I will expand this later).
#'   - Adds residual noise, and scales those residuals so heritability is
#'       known exactly (the script loops over h2 values, hard-coded as 0.2,
#'       0.4, 0.6 and 0.8).
#'   - Saves a phenotype data file that GEMMA can read.
#'
#' It is assumed that input matrices contain the same sets of SNPs for all
#' three cohorts. Effect sizes are reused for the three cohorts so they are
#' directly comparable.
#'
#' It is also assumed that the SNP matrices contain N SNPs, to generate
#' k replicate datasets of n SNPs each. The script will return an error
#' if N/n is not a round number.

library("optparse")
suppressPackageStartupMessages(library("tidyverse"))

# === Define functions ===

#' Import genotype data
#'
#' Import a matrix of SNPs, transpose it, and standardise variances at each locus
#' @param file_path str Path to a tab-separated text file with a row for each
#' locus and a column for each sample.
import_and_transpose_matrix <- function(file_path){
  geno_matrix <- read_tsv(file_path, show_col_types = FALSE)
  # geno_matrix has a row for each marker and a column for each line.
  # Transpose it for a more logical tidy table.
  geno_matrix <- geno_matrix %>%
    pivot_longer(-snp, names_to = "line", values_to = "value") %>%
    pivot_wider(names_from = snp, values_from = value)
  # Standardise columns
  stand_geno_matrix <- apply(geno_matrix[,-1], 2, standardise_genetic_variance) %>%
    as_tibble()

  return(list(
    line_names = geno_matrix$line,
    geno = stand_geno_matrix)
  )
}

# Check for columns with no variance or lots of missing data.
check_input_matrix <- function(geno_matrix){
  for(col in 2:ncol(geno_matrix)){
    column_var <- var(geno_matrix[,col], na.rm = TRUE)
    if(column_var == 0){
      cat(names(geno_matrix)[col], "has zero variance.\n")
    }
    frequency_NA <- mean(is.na(geno_matrix[,col]))
    if(frequency_NA > 0.1){
      cat(names(geno_matrix)[col], "has more than 10% NAs.\n")
    }
  }
}

# Function to standardise variance at each locus so effects are orthogonal.
standardise_genetic_variance <- function(allele_counts){
  # Impute missing data
  # model.matrix will silently drop anything with NAs
  # That will be really problematic when there are many loci.
  # To get around this, fill in missing data with the mean genotype value.
  # This is what GEMMA will do anyway later on.
  allele_counts[is.na(allele_counts)] <- mean(allele_counts, na.rm = TRUE)
  # SNPs were picked to have non-zero MAC *in the parents*, but sometimes
  # genetic variants have no variance.
  # This is especially true for the crosses, because there can be some fixation,
  # even if the parents had variance at the same locus.
  # In this case, I still want to keep them so I can properly compare between
  # cohorts, so set the variance to zero.
  genetic_variance <- var(allele_counts, na.rm=TRUE)
  if(genetic_variance == 0){
    return(allele_counts * 0)
  } else {
    # If there is finite variance, perform the standardisation.
    p <- mean(allele_counts/2, na.rm=TRUE)
    standardised_values <- (allele_counts-2*p) / sqrt(2*p*(1-p))

    return(standardised_values)
  }
}



#' Get indices to split up a genotype matrix
#'
#' The genotype matrix contains m*k=N columns, and I want to split it into
#' m sub-matrices each of size k.
#' This function returns a list of length m, each with k integer indices.
#' @param geno_matrix Output of `import_and_transpose_matrix`
#' @param loci_per_model Number of loci (k) per subset.
get_column_indices_to_subset_input <- function(geno_matrix, loci_per_model){
  total_loci <-ncol(geno_matrix)
  # If m*k=N, create a list of indices
  if(total_loci %% loci_per_model == 0){
    n_subsets <- total_loci / loci_per_model
    col_ix <- split(1:total_loci, cut(seq_along(1:total_loci), n_subsets, labels = FALSE))
    return(col_ix)
  } else {
    # If m*k =/= N, return an error.
    cat("The total number of loci in the file is not a multiple of the number of loci for each model.")
    stop()
  }
}


# Simulate genetic values and add residual noise.
# This ensures heritability is *exactly* what you expect it to be.
simulate_phenotypes <- function(model_matrix, effect_sizes, heritability, line_names){
  # Simulate expected values.
  genetic_values <- mm %*% effect_sizes
  genetic_variance <- as.numeric(var(genetic_values))
  # Normally distributed residuals for each line
  residuals <- rnorm(length(genetic_values))
  # Ensure the finite sample really has variance = st deviation = 1
  residuals <- residuals / sd(residuals, na.rm=TRUE)
  # Scale the residuals in proportion to heritability
  scaling_factor <- sqrt(genetic_variance * (1/heritability -1))
  residuals <- residuals  * scaling_factor
  # Create a tibble in the form GEMMA can use.
  phenotypes <- tibble(
    FID=0, # column of zeroes
    line = line_names, # Column of sample names
    phenotype = as.numeric(genetic_values + residuals) # Column of phenotypes
  )

  return(phenotypes)
}

# === Parse arguments ===

parser <- OptionParser()
parser <- add_option(parser, "--parents",
                     help="Genotype matrix file for the parents.")
parser <- add_option(parser, "--cohort1",
                     help="Genotype matrix file for F8 cohort 1.")
parser <- add_option(parser, "--cohort2",
                     help="Genotype matrix file for F8 cohort 2.")
parser <- add_option(parser, "--nloci",
                     type = "integer",
                     help="Number of loci in each model.")
parser <- add_option(parser, "--output",
                     help="Output prefix for three phenotype files.")
opt <- parse_args(parser)

# Ensure the directory ends in a /
outdir <- paste0(fs::path_norm(opt$output),"/")

cat("Starting simulations using", opt$nloci, "loci for genotype matrices:\n")
cat(opt$parents, "\n")
cat(opt$cohort1, "\n")
cat(opt$cohort2, "\n")
cat("\nOutput files will be saved to", opt$output)
cat("\n\n")

# === Main ===

input_files<- list(
  parents = opt$parents,
  cohort1 = opt$cohort1,
  cohort2 = opt$cohort2
)

# Import the data as a list with an element for each cohort.
genotype_matrix_list <- lapply(input_files, import_and_transpose_matrix)

# Ensure the matrices have the same number of loci (columns)
for(i in 1:(length(genotype_matrix_list)-1)){
  for(j in 2:length(genotype_matrix_list)){
    if(j>i){
      dim_i <- ncol(genotype_matrix_list[[i]]$geno)
      dim_j <- ncol(genotype_matrix_list[[j]]$geno)
      if(dim_i != dim_j){
        cat("Error: genotype matrices do not have the same numbers of loci.")
        stop()
      }
    }
  }
}
# List of vectors of indices to subset the genotype matrix
col_ix <- get_column_indices_to_subset_input(
  genotype_matrix_list$parents$geno,
  loci_per_model = opt$nloci
)

# For each set of markers, draw effect sizes, simulate phenotypes and
# save a phenotype file to disk.
# One set of effect sizes are drawn and applied to all three cohorts.
cat("\nStarting simulations for replicate ")
for(i in 1:length(col_ix)){
  cat(i, "")
  # For this simulation, effect sizes don't really matter, because
  # there is only one locus but I will include them for completeness,
  # and in case I change this later.
  effect_sizes <- rnorm(opt$nloci)

  # Use the effect sizes to simulate phenotypes for the three cohorts.
  for(cohort in names(genotype_matrix_list)){
    cohort_genotypes <- genotype_matrix_list[[cohort]]$geno[,col_ix[[i]]]
    # Model matrix for the regression
    mm <- model.matrix(~ 0 + (.), data=cohort_genotypes)

    # Loop over heritability values
    for(heritability in c(0.2, 0.4, 0.6, 0.8)){
      # Create a table of phenotypes in a format GEMMA can read.
      phenotype_table <- simulate_phenotypes(
        model_matrix = mm,
        effect_sizes = effect_sizes,
        heritability = heritability,
        line_names = genotype_matrix_list[[cohort]]$line_names
      )

      # A file path containing cohort, replicate, MAF range and heritability.
      rep_id <- sprintf("%03d", i)
      maf_range <- str_extract(input_files[[cohort]], "0.[0-3]-0.[0-4]")
      outfile_path <- paste0(
        outdir, cohort, "_", rep_id, "_", maf_range, "_", heritability, ".tsv"
      )
      # Write to disk.
      write_tsv(
        phenotype_table,
        file = outfile_path,
        col_names = FALSE
      )
    }
  }
}
cat("\nFinished.\n\n\n")

