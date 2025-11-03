#' Get GPS coordinates for each parental line.
#'
#' Most are taken from an updated file from the 1001 genomes project
#' The remaining 18 accessions are from Anastasio et al 2011.

library(tidyverse)

# A list of names for parental accessions.
parental_names <- read_csv(
  "03_processing/09_impute_haplotypes/output/parental_line_names.txt",
  col_names = c("id"),
  col_types = 'c'
)

# Metadata for the 1135 accessions in the 1001 genomes paper, with corrections
# by Tal Dahan.
metadata_1135 <- read_csv(
  "01_data/09_parental_metadata/dahan_updated_metadata.csv",
  show_col_types = FALSE
) %>%
  mutate(
    id = as.character(id)
  ) %>%
  select(
    id, latitude, longitude, group
  )
# Metadata for 5966 accessions from Anastasio et al 2011
metadata_5966 <- read_csv(
  "01_data/09_parental_metadata/anastasio_et_al_2011.csv",
  show_col_types = FALSE
) %>%
  rename(
    id = `AccessionNo.`,
    latitude = Latitude,
    longitude = Longitude
  ) %>%
  mutate(
    id = as.character(id),
    group = NA
    ) %>%
  select(
    id, latitude, longitude, group
  )



# Get coordinates for most of the parents
parental_coords <- parental_names %>%
  left_join(metadata_1135, by = "id")

# 18 lines are not in the 1135
missing_ids <- parental_coords %>%
  filter(is.na(latitude)) %>%
  pull(id)

# Bind data on accessions from the two metadatasets
parental_coords <- rbind(
  parental_coords %>%
    filter(! id %in% missing_ids),
  metadata_5966 %>%
    filter(id %in% missing_ids)
)


parental_coords <- parental_coords %>%
  mutate(
    latitude = case_when(
      id == "7519" ~ 56.14,
      id == "5829" ~ 55.382301,
      id == "8423" ~ 56.113287,
      .default = latitude
    ),
    longitude = case_when(
      id == "7519" ~ 15.78,
      id == "5829" ~ 14.053164,
      id == "8423" ~ 13.729755,
      .default = longitude
    )
  )


parental_coords %>%
  write_csv(
    "03_processing/06_process_phenotypes/output/parental_GPS_coordinates.csv"
    )

