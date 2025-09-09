#' Join GEMMA results files for resampled datasets.
#'
#' Create text files summarising pvalues and effect sizes from each resampled
#' dataset. Columns indicate chr, position, then many columns giving test
#' statistics from each dataset. This is done separately for analyses with and
#' without a K matrix correction.
#'
#' Tom Ellis, 23rd June 2025

suppressPackageStartupMessages(
  library(tidyverse)
)

# Vector of paths to GEMMA results files
dataset_paths <- Sys.glob(
  "/scratch-cbe/users/thomas.ellis/crosses/19_ft12_resampled/ft12_resample_*"
)

# Import and merge results files.
# This is done separately for files with and without a K matrix.
for(k in c("no_K", "with_K")){
  cat("Processing files", k, "\n")
  # Initializes the progress bar
  pb <- txtProgressBar(0, length(dataset_paths))

  for(i in 1:length(dataset_paths)){
    setTxtProgressBar(pb, i)

    # Path to GEMMA results file
    path_prefix = "/scratch-cbe/users/thomas.ellis/crosses/19_ft12_resampled/ft12_resample_"
    assoc_file <- paste0(path_prefix, i, "/",k,"/phenotypes.assoc.txt" )
    # Import data and calculate -log10 pvalues
    pvals_i <- read_delim(
      assoc_file,
      delim = "\t",
      show_col_types = FALSE
      ) %>%
      mutate(
        log10p = -log10(p_lrt)
      )
    # Copy the dataset, for saving effect sizes and allele frequencies.
    betas_i <- afreq_i <- pvals_i

    # Select only chr, position, and either pvalues, effect sizes or allele frequency
    pvals_i <- pvals_i %>% select(chr, ps, log10p)
    betas_i <- betas_i %>% select(chr, ps, beta)
    afreq_i <- afreq_i %>% select(chr, ps, af)
    # Change the last column to match the name of the dataset
    names(pvals_i)[3] <- paste0("log10p_", i)
    names(betas_i)[3] <- paste0("beta_", i)
    names(afreq_i)[3] <- paste0("afreq_", i)

    # If this is the first dataset, just copy the results
    if(i == 1){
      merge_pvals <- pvals_i
      merge_betas <- betas_i
      merge_afreq <- afreq_i
      # Create a column for getting the metastatistics
      merge_pvals <- merge_pvals %>% mutate(sum_log_p = log10p_1)
      merge_betas <- merge_betas %>% mutate(sum_beta  = beta_1)
      merge_afreq <- merge_afreq %>% mutate(sum_afreq = afreq_1)
    }
    # If this is not the first data set, perform an outer join to
    # previous results.
    # Add the result to the running sum.
    if(i > 1){
      merge_pvals <- merge_pvals %>%
        full_join(pvals_i, by = c('chr', 'ps')) %>%
        mutate(
          across(where(is.numeric), ~replace_na(.x, 0)),
          sum_log_p = sum_log_p + eval(parse(text=paste0("log10p_", i)))
        )
      merge_betas <- merge_betas %>%
        full_join(betas_i, by = c('chr', 'ps')) %>%
        mutate(
          across(where(is.numeric), ~replace_na(.x, 0)),
          sum_beta = sum_beta + eval(parse(text=paste0("beta_", i)))
        )
      merge_afreq <- merge_afreq %>%
        full_join(afreq_i, by = c('chr', 'ps')) %>%
        mutate(
          across(where(is.numeric), ~replace_na(.x, 0)),
          sum_afreq = sum_afreq + eval(parse(text=paste0("afreq_", i)))
        )
    }
  }
  close(pb) # Close the connection

  # Round to 3 decimal places to save disk space.
  merge_pvals <- merge_pvals %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
  merge_betas <- merge_betas %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
  merge_afreq <- merge_afreq %>%
    mutate(across(where(is.numeric), ~ round(., 3)))

  # Write objects to disk.
  # Files with stats for each dataset.
  merge_pvals %>%
    select(-sum_log_p) %>%
    write_csv(
      paste0("05_results/19_ft12_resampled/output/ft12_resampled_pvals_", k, ".csv")
    )
  merge_betas %>%
    select(-sum_beta) %>%
    write_csv(
      paste0("05_results/19_ft12_resampled/output/ft12_resampled_betas_", k, ".csv")
    )
  merge_afreq %>%
    select(-sum_afreq) %>%
    write_csv(
      paste0("05_results/19_ft12_resampled/output/ft12_resampled_afreq_", k, ".csv")
    )


  # Files with metastatistics
  merge_pvals %>%
    select(chr, ps, sum_log_p) %>%
    write_csv(
      paste0(
        "05_results/19_ft12_resampled/output/ft12_resampled_sum_pvals_", k, ".csv")
    )
  merge_betas %>%
    mutate(
      mean_beta = sum_beta / length(dataset_paths)
    ) %>%
    select(chr, ps, mean_beta) %>%
    write_csv(
      paste0("05_results/19_ft12_resampled/output/ft12_resampled_mean_betas_", k, ".csv")
    )
  merge_afreq %>%
    mutate(
      mean_afreq = sum_afreq / length(dataset_paths)
    ) %>%
    select(chr, ps, mean_afreq) %>%
    write_csv(
      paste0("05_results/19_ft12_resampled/output/ft12_resampled_mean_afreq_", k, ".csv")
    )
}

