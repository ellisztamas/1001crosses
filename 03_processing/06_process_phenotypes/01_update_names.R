library(tidyverse)

# ibdpainting results on the resequenced dataset
reseq <- read_tsv(
  "03_processing/04_resequencing/output/ibdpainting_results_resequenced.tsv"
  ) %>%
  dplyr::select(name, new_name) %>%
  dplyr::filter(! name %in% c("1137", '1074')) %>%  # Exclude parental lines
  dplyr::filter(! grepl("blank", name))  # Exclude blank wells

# ibdpainting results on the low-coverage dataset
low_coverage <- read_csv(
  "03_processing/03_validate_genotypes/output/ibdpainting_results.csv"
  ) %>%
  mutate(
    name     = gsub(" F8 ", "_", name),
    new_name = gsub(" F8 ", "_", new_name) # Change new names
  ) %>%
  dplyr::filter(! name %in% reseq$name)  %>%  # Exclude lines that were resequenced
  dplyr::filter(! grepl("Aa", name)) %>%  # Exclude A. arenosa
  # There are two rows with name==6012x997_rep1
  # One corrects to 6021x9407_rep1, the other does not.
  # I think someone messed up the sequence sheet, but the phenotype data are fine
  dplyr::filter( new_name != "6021x9407_rep1")  %>%
  dplyr::select(name, new_name)


naming_file <- rbind(reseq, low_coverage) %>%
  arrange(new_name) %>%
  dplyr::rename(
    original_name = name,
    corrected_name= new_name
  )
