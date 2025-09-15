library(dplyr)

#' Calculate emission probabilities.
#'
#' Function to calculate emission probabilities that an F8 individual
#'
#' @param vcf_file VCF file in vcfR format. Import this with `read.vcfR()` from
#' vcfR.
#' @param parent1 String giving the name of the first parent.
#' This should be one of the samples in `vcf_file`.
#' @param parent2 String giving the name of the second parent.
#' This should be one of the samples in `vcf_file`.
#' @param progeny String giving the name of the progeny individual.
#' This should be one of the samples in `vcf_file`.
#' @param err Probability that a read aligns in the wrong place.
#'
#' @return Tibble giving chromosome, base position, distance from the previous
#' marker, and emission probabilities that the ancestry at each marker is from
#' parent1, heterozygous or parent2.
#'
#' @author Tom Ellis
#'
#' Calculate emission probabilities for each locus that the progeny genotype
emission_probabilities <- function(vcf_file, parent1, parent2, progeny, err){

  chrom <- getCHROM(vcf_file)
  pos <- getPOS(vcf_file)


  cat("Extracting (hard) genotype calls in the parents.\n\n")
  # Matrix of (hard) genotype calls
  gt <- extract.gt(vcf_file, element = "GT")
  # If there are any phased genotypes, remove phasing information
  # for consistent comparisons between samples
  gt[,parent1] <- gsub("\\|", "/", gt[,parent1])
  gt[,parent2] <- gsub("\\|", "/", gt[,parent2])
  gt[,progeny] <- gsub("\\|", "/", gt[,progeny])
  # Make heterozygous calls consistent
  gt[,parent1] <- gsub("1///0", "0/1", gt[,parent1])
  gt[,parent2] <- gsub("1\\/0", "0/1", gt[,parent2])

  cat("Extact offspring genotype-call probabilities from VCF field 'AD'.\n\n")
  # Extract genotype likelihoods as strings (eg. "13,0")
  ad_string <- extract.gt(vcf_file, element = "AD")
  # Convert missing data to zeroes, for better or worse
  ad_string <- gsub("\\.", "0", ad_string[,progeny])
  # List of pairs
  ad_split <- strsplit(ad_string, ",")
  # Convert the list of character vectors to a list of numeric vectors
  ad_numeric_list <- lapply(ad_split, as.numeric)
  # 3. Extract the first element (REF depth) from each sublist
  ref_depths <- sapply(ad_numeric_list, `[`, 1)
  # 4. Extract the second element (ALT depth) from each sublist
  alt_depths <- sapply(ad_numeric_list, `[`, 2)

  # Probabilities that markers are each of three genotypes.
  call_probs <- tibble(
    hom0 = (1-err)^ref_depths * err^(alt_depths),
    het  = 0.5^ref_depths * 0.5^alt_depths,
    hom1 = (1-err)^alt_depths * err^(ref_depths)
  )


  # Calculate emission probabilties for each locus.
  # This returns a list of vectors with three elements giving the probabilties
  # that F8 data are homozygous parent 1, heterozygous, or homozygous parent 2.
  cat("Calculating emission probabilties at each locus.\n\n")
  emission_prob_list <- vector('list', nrow(gt))
  # Initializes the progress bar=
  pb <- txtProgressBar(min = 0, max = nrow(gt))
  for(i in 1:nrow(gt)){
    setTxtProgressBar(pb, i)

    # Parental genotype calls and offspring genotype probs for a single row.
    pg1 <- gt[i,parent1]
    pg2 <- gt[i,parent2]
    ofg <- as.numeric(call_probs[i,])
    # Ignore
    if( is.na(pg1) | is.na(pg2) ) next
    # Ignore cases where genotypes are identical
    if(pg1 == "0/0" & pg2 == "0/0") next
    if(pg1 == "1/1" & pg2 == "1/1") next
    # Opposing homozygotes
    if(pg1 == "0/0" & pg2 == "1/1") p <- ofg
    if(pg1 == "1/1" & pg2 == "0/0") p <- ofg[3:1]
    # When one parent is heterozygous, weight towards the other parent
    if(pg1 == "0/0" & pg2 == "0/1") p <- ofg      * c(2,1,1)
    if(pg1 == "1/1" & pg2 == "0/1") p <- ofg[3:1] * c(2,1,1)
    if(pg1 == "0/1" & pg2 == "0/0") p <- ofg[3:1] * c(1,1,2)
    if(pg1 == "0/1" & pg2 == "1/1") p <- ofg      * c(1,1,2)
    # When both parents are heterozygous, any thing is possible
    if(pg1 == "0/1" & pg2 == "0/1") p <- ofg + ofg[3:1]
    # Normalise to sum to one
    p <- p / sum(p)
    p <- c(i, p)
    names(p) <- c('i', 'parent1', 'het', 'parent2')
    emission_prob_list[[i]] <- c(p)
  }
  close(pb)

  cat("Merging SNP positions.\n\n")
  # Add marker positions, remove uniformative markers
  emission_data <- do.call('rbind', emission_prob_list) %>%
    as_tibble() %>%
    mutate(
      chr = chrom[.$i],
      pos = pos[.$i]
    ) %>%
    filter(
      complete.cases(.)
    )

  # Add distances between markers (separately for each chromosome)
  chr_labels <- unique(emission_data$chr)
  dist <- vector('list', length(chr_labels))
  for(chr in chr_labels){
    chr_ix <- emission_data$chr == chr
    chr_pos <- emission_data$pos[chr_ix]
    dist[[chr]] <- c(0, chr_pos[2:length(chr_pos)] - chr_pos[1:(length(chr_pos)-1)])
  }
  emission_data$dist <- do.call('c', dist)

  # Rearrange columns for a more logical order.
  emission_data <- emission_data %>%
    select(chr, pos, dist, parent1, het, parent2)


  return(emission_data)
}
