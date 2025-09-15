library(dplyr)


#' Calculate most-likely hidden states
#'
#' Use the Viterbi algorithm to calculate the most likely hidden state at each
#' marker.
#'
#' @param emission_data Output of `emission_probabilities`.
#' @param remove_singletons Bool. If TRUE, haplotypes consisting of only a
#' single marker are removed from the BED file output.
#'
#' @return A list of two tibbles.
#' `markers` returns the input tibble with an additional column showing
#' most-likely hidden states.
#' `haplotypes` returns a tibble summaring runs of hidden states as haplpotyes.
#' Columns indicate chromosome, start and end points of each haplotype, the
#' ancestry at each haplotype, and the number of markers included.
#'
#' @author Tom Ellis
viterbi_hmm <- function(emission_data, remove_singletons=TRUE) {
  # Extract emission probabilities and distances
  emission_probs <- emission_data[, c("parent1", "het", "parent2")]
  chr_labels <- unique(emission_data$chr)

  # Define states
  states <- c("parent1", "het", "parent2")
  n_states <- length(states)
  n_obs <- nrow(emission_data)
  cat("Calculating most likely path through", n_obs, "steps with", n_states, "hidden states.\n")

  # Calculate transition probabilities
  # Distance in bp -> distance in megabases
  # 4cM per megabase
  # 100 cM per Morgan
  transition_rate <- emission_data$dist*1e-6 * 4 * 1e-2

  # Initialize transition probability matrices for each step
  transition_probs <- array(0, dim = c(n_obs - 1, n_states, n_states))

  cat("Calculating transition probabilities at each step.\n")
  for (t in 1:(n_obs - 1)) {
    # Off-diagonal elements (between different states)
    for (i in 1:n_states) {
      for (j in 1:n_states) {
        if (i != j) {
          transition_probs[t, i, j] <- transition_rate[t] / 2  # Split equally between other states
        }
      }
    }
    # Diagonal elements (staying in same state)
    for (i in 1:n_states) {
      transition_probs[t, i, i] <- 1 - transition_rate[t]
    }
  }

  cat("Running the Forward-Backward algorithm to determine the most likely states.\n")
  # Initialize Viterbi tables
  viterbi_prob <- matrix(0, nrow = n_obs, ncol = n_states)
  viterbi_path <- matrix(0, nrow = n_obs, ncol = n_states)

  # Initialize first observation (assuming uniform prior)
  viterbi_prob[1, ] <- log(1/n_states) + log(as.numeric(emission_probs[1, ]))

  # Forward pass
  for (t in 2:n_obs) {
    for (j in 1:n_states) {
      # Calculate probabilities from all previous states
      trans_probs <- viterbi_prob[t-1, ] + log(transition_probs[t-1, , j])

      # Find maximum and store
      max_idx <- which.max(trans_probs)
      viterbi_prob[t, j] <- trans_probs[max_idx] + log(as.numeric(emission_probs[t, j]))
      viterbi_path[t, j] <- max_idx
    }
  }

  # Backward pass - find most likely path
  path <- numeric(n_obs)
  path[n_obs] <- which.max(viterbi_prob[n_obs, ])

  for (t in (n_obs-1):1) {
    path[t] <- viterbi_path[t+1, path[t+1]]
  }

  # Convert to state names
  emission_data$ancestry <- states[path]

  # Remove
  if(remove_singletons){
    cat('Purging haplotypes that are visited for only a single marker.\n')
    corrected_state <- vector('list', length(chr_labels))
    names(corrected_state) <- chr_labels

    for(chr in chr_labels){
      chr_state <- emission_data$ancestry[emission_data$chr == chr]
      n <- length(chr_state)

      singleton_ix <- which(
        c(chr_state[1:(n-1)]    != chr_state[2:n], FALSE) &
          c(FALSE, chr_state[2:n] != chr_state[1:(n-1)])
      )

      chr_state[singleton_ix] <- chr_state[singleton_ix-1]
      corrected_state[[chr]] <- chr_state
    }
    emission_data$ancestry <- do.call('c', corrected_state)
  }

  cat("Calculating the start and end points of each haplotype.\n")
  haplotypes_list  <- vector('list', length(chr_labels))
  names(haplotypes_list) <- chr_labels

  # Known lengths of the TAIR10 assembly
  chr_lengths <- list(
    Chr1 = 30427671,
    Chr2 = 19698289,
    Chr3 = 23459830,
    Chr4 = 18585056,
    Chr5 = 26975502
  )

  for(chr in chr_labels){
    chr_states <- emission_data[emission_data$chr == chr,]
    rle_states <- rle(chr_states$ancestry)
    # Calculate start and end indices for each run
    ends <- cumsum(rle_states$lengths)
    starts <- ends - rle_states$lengths + 1

    # Get start and end positions from the original data
    chr_haplotypes <- tibble(
      chr = chr,
      start_position = chr_states$pos[starts],
      end_position = chr_states$pos[ends],
      state = rle_states$values,
      n_markers = rle_states$lengths
    )
    # Fill in the start and ends of the chromosomes.
    chr_haplotypes$start_position[1] <- 1
    chr_haplotypes$end_position[nrow(chr_haplotypes)] <- chr_lengths[[chr]]

    haplotypes_list[[chr]] <- chr_haplotypes
  }
  haplotypes <- do.call('rbind', haplotypes_list)


  return(list(
    markers = emission_data,
    haplotypes = haplotypes
  ))
}

