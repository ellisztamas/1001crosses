#' Format flowering time data collected by Almudena
#'
#' Alumdena Molla Morales recorded flowering time data for two weeks in May and
#' July 2024 while Tom Ellis was away. She used a different app to scan bar codes
#' and saved the data into separate files for each data.
#'
#' This script parses those files and creates a dataframe in the same format as
#' the data from MementoDB
#'
#' Input:
#'    - 24 files listing plants that flowered on a single day
#' Output:
#'    - Dataframe with a row for each plant and the following columns:
#'        key : unique key for each plant.
#'        qr : Text from the QR code on each plant listing unique ID, genotype and position
#'        date : Data of flowering
#'        user : All "Almudena Molla Morales"
#'        flowers_tomorrow : All False for these data
#'        flowered_yesterday : All False for these data
#'        comment: All empty for these data

library('tidyverse')

files_by_day <- c(
  Sys.glob("01_data/05_phenotype_expt/flowering_time_raw_data/*2024053134631.csv"),
  Sys.glob("01_data/05_phenotype_expt/flowering_time_raw_data/*sun*613.csv"),
  Sys.glob("01_data/05_phenotype_expt/flowering_time_raw_data/*655.csv"),
  Sys.glob("01_data/05_phenotype_expt/flowering_time_raw_data/*502.csv"),
  Sys.glob("01_data/05_phenotype_expt/flowering_time_raw_data/*028.csv")
)
ft10_AMM <- vector('list', length(files_by_day))

for(d in 1:length(files_by_day)){
  today <- files_by_day[d]

  date <- paste(
    substr(basename(today), 7, 8),
    substr(basename(today), 5, 6),
    substr(basename(today), 1, 4),
    sep="/"
  )

  qr <- read_csv(today, show_col_types = FALSE) %>%
    filter(`#` != "Total") %>%
    pull(Barcode)

  ft10_AMM[[d]] <- tibble(
    key = NA,
    qr = qr,
    date = date,
    user = "Almudena Molla Morales",
    flowers_tomorrow = FALSE,
    flowered_yesterday = FALSE,
    comment = NA,
    did_not_germinate = FALSE,
    did_not_flower = FALSE
  )
}

ft10_AMM <- do.call('rbind', ft10_AMM)
