#' Fill in gaps in .bed files
#'
#' Previous scripts created .bed files giving identifiable start and end
#' points of parental haplotypes in each F8, as well as some heterozygous and
#' private haplotypes. However, these haplotypes are not all contiguous.
#'
#' This script takes each .bed file and adds additional lines to fill in those
#' gaps. For each chromosome separately it identifies
#'
#' Tom Ellis 21st August 2025

library(tidyverse)

# Vector of input .bed files to be processed.
validated_bed_files <- Sys.glob(
  "03_processing/09_impute_haplotypes/output/08_check_breakpoints/*bed"
)
# There are ten bed files to be skipped.
bed_files_to_skip <- c(
  '9381x6198_rep1', # This was repeated as 9381x6198
  '9381x6198_rep2', # This was repeated as 9381x6198
  '5829x9336_rep1', # This, and those below, were identified as having incorrect parentage
  '5829x9336_rep2',
  '6030x992_rep2',
  '8283x9471_rep1',
  '8283x9471_rep2',
  '9380x1063_rep1',
  '9399x6097_rep1',
  '9399x6097_rep2'
)
validated_bed_files <- grep(
  pattern = paste(bed_files_to_skip, collapse = "|"),
  x = validated_bed_files,
  value = TRUE,
  invert = TRUE
  )

# Directory to store new .bed files with additional rows.
outdir <- "03_processing/09_impute_haplotypes/output/09_fill_haplotype_gaps"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# A log file to catch any issues.
logfile <- paste0(outdir, "/", "logfile.txt")
cat(paste0(date(), "\n"), file=logfile, append=FALSE)

for(file_path in validated_bed_files){
  # Load the .bed file
  sample_name <- basename(file_path)
  input_bed_file <- read_tsv(file_path, show_col_types = FALSE)
  chr_names <- unique(input_bed_file$chr)
  cat("Processing file", sample_name, "\n")
  # Empty list to store new bed data for each chromosome.
  output_bed_file <- vector('list', length(chr_names))
  names(output_bed_file) <- chr_names

  # Create new bed files for each chromosome.
  for(chr_label in chr_names){

    if(!any(input_bed_file$chr == chr_label)){
      cat(
        chr_label, "is missing entirely from.", sample_name, "\n",
        file = logfile, append=TRUE
      )
    }
    # Extract a single chromosome.
    chr_bed <- input_bed_file %>%
      filter(chr == chr_label) %>%
      arrange(start_position)


    # Check the bed file contains position 1.
    if(! 1 %in% chr_bed$start_position){
      cat(
        chr_label, "is missing postiion 1 in", sample_name, "\n",
        file = logfile, append=TRUE
      )
    }

    # Check the bed file extends to the end of the chromosome.
    expected_chr_lengths <- list(
      Chr1 = 30427671,
      Chr2 = 19698289,
      Chr3 = 23459830,
      Chr4 = 18585056,
      Chr5 = 26975502
    )
    if(! expected_chr_lengths[[chr_label]] %in% chr_bed$end_position){
      cat(
        chr_label, "is missing a final postiion in", sample_name, "\n",
        file = logfile, append=TRUE
      )
    }

    # Check for consecutive matching states.
    if(any(rle(chr_bed$state)$lengths > 1)){
      cat(
        "One or more states on", chr_label, "in",
        sample_name, "contains consecutive matching states.\n",
        file = logfile, append=TRUE
      )
    }
    # Check for any haplotypes that overlap
    any_overlaps <- any(
      chr_bed$start_position[-1] < chr_bed$end_position[1:(nrow(chr_bed)-1)]
    )
    if(any_overlaps){
      cat(
        "Two or more haplotypes on", chr_label, "in",
        sample_name, "overlap.\n",
        file = logfile, append=TRUE
      )}

    # Add additional rows between haplotypes.
    if(nrow(chr_bed) > 1){
      for(r in 1:(nrow(chr_bed)-1)){
        distance_between_haplpotypes <- chr_bed$start_position[r+1] - chr_bed$end_position[r]
        if(distance_between_haplpotypes > 1){
          new_row <- tibble(
            chr = chr_label,
            start_position = chr_bed$end_position[r] + 1,
            end_position   = chr_bed$start_position[r+1] -1,
            state = "private",
            n_markers = NA
          )
          chr_bed <- rbind(chr_bed, new_row)
        }
      }
    }
    # Send data for this chromosome to output_bed_file
    output_bed_file[[chr_label]] <- chr_bed %>%
      arrange(start_position)
  }
  # Bind the data for each chromosome and write to disk.
  do.call('rbind', output_bed_file) %>%
    write_tsv(paste0(outdir, "/", sample_name))
}
