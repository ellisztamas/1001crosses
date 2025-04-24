library('tidyverse')

ft10_TJE <- read_csv("01_data/05_phenotype_expt/flowering_time_raw_data/241018_ft10.csv")

source("01_data/05_phenotype_expt/format_ft10_from_AMM.R")

# For 36 it was not possible to read the QR code, so only the integer key was
# recorded.
# Fix these manually.
flowering_time <- rbind(ft10_AMM, ft10_TJE) %>%
  mutate(
    qr = case_match(
      qr,
      '1000' ~ '1000_1.1.A6_9323x8369 rep1_F8',
      '1012' ~ '1012_1.8.B6_9336x6220 rep1_F8',
      '1024' ~ '1024_1.1.F3_9352x8422 rep1_F8',
      '1054' ~ "1054_1.3.E2_9380x1063 rep1_F8",
      "1188" ~ "1188_3.37.G3_9433x6030 rep1_F8",
      "1207" ~ "1207_1.11.C1_9450x9356 rep2_F8",
      "1219" ~ "1219_1.9.C4_9452x6114 rep2_F8",
      "124"  ~ "124_1.10.C3_5835x9058 rep2_F8",
      "1258" ~ "1258_1.3.G6_9481x9399 rep1_F8",
      "1440" ~ "1440_3.38.E2_6086_parent",
      "163"  ~ "163_1.3.E5_6017x9392 rep2_F8",
      "187"  ~ "187_1.7.B4_6024x1367 rep1_F8",
      "190"  ~ "190_1.1.H4_6024x1367 rep2_F8",
      "214"  ~ "214_1.2.C2_6034x6231 rep2_F8",
      "253"  ~ "253_1.9.F2_6043x9453 rep1_F8",
      "280"  ~ "280_1.6.D5_6070x6019 rep2_F8",
      "313"  ~ "313_1.7.F5_6086x8222 rep2_F8",
      "322"  ~ "322_1.2.B1_6090x9416 rep1_F8",
      "409"  ~ "409_1.3.C5_6115x1585 rep2_F8",
      "421" ~ "421_1.4.H2_6124x6184 rep1_F8",
      "430" ~ "430_1.1.E5_6125x6099 rep2_F8",
      "451" ~ "451_1.10.H4_6131x6113 rep1_F8",
      "475" ~ "475_1.8.G5_6134x6013 rep1_F8",
      "541" ~ "541_1.12.A1_6154x6221 rep1_F8",
      "589" ~ "589_1.8.D3_6177x6023 rep1_F8",
      "631" ~ "631_1.2.A3_6193x9394 rep1_F8",
      "643" ~ "643_1.6.D1_6195x1074 rep1_F8",
      "664" ~ "664_1.2.C3_6201x9399 rep2_F8",
      "724" ~ "724_1.6.H5_6237x6090 rep1_F8",
      "799" ~ "799_1.4.A2_6913x6149 rep2_F8",
      "820" ~ "820_1.6.G2_6974x6102 rep1_F8",
      "853" ~ "853_1.4.A4_8227x6216 rep2_F8",
      "904" ~ "904_1.12.C6_8258x6115 rep1_F8",
      "922" ~ "922_1.2.A1_8307x9481 rep1_F8",
      "943" ~ "943_1.8.F6_8369x6132 rep2_F8",
      "991" ~ "991_1.5.E3_9058x6112 rep2_F8",
      .default= qr
    )
  )


# Three plants flowered on the 1st June and need the date correcting.
flowering_time <-  flowering_time %>%
  mutate(
    date = ifelse( grepl("1st June",    comment), dmy("01/06/2024"), date)
  )


flowering_time <- flowering_time %>%
  filter(
    qr != "pot with no tag from tray 9" # Remove row with a string for a code.
    # key == 173 # This plant is given twice. I think this entry is wrong, bc the days to flower is much bigger than the other reps
  ) %>%
  # Three rows have one underscore too many in their code.
  mutate(qr = str_replace_all(qr, "_rep", " rep")  ) %>%
  # Split the code
  separate(qr, into=c("id", "position", "genotype", "generation"), sep = "_") %>%
  separate(position, into=c("replicate", 'tray', 'position'), remove=TRUE) %>%
  mutate(
    tray = as.integer(tray),
    # Column to separate cohort groups of F9s and the parents
    cohort = case_when(
      str_detect(genotype, "rep1") ~ "rep1",
      str_detect(genotype, "rep2") ~ "rep2",
      generation == "parent" ~ "parents"
    ),
    # date the plant flowered
    date = dmy(date),
    # Date the plant was sown (or moved to the growth room)
    start_date = dmy(
      case_when(
        tray %in% 1:14 ~ "27/03/2024",
        tray %in% 15:20 ~ "16/04/2024",
        tray %in% 21:28 ~ "16/04/2024",
        tray %in% 29:34 ~ "22/04/2024",
        tray %in% 35:42 ~ "23/04/2024"
      )
    ),
    # Number of days to flower.
    days_to_flower = as.integer(
      difftime(date, start_date, units = "days")
    ),
    # Adjust for being out by a day.
    days_to_flower = days_to_flower - flowered_yesterday + flowers_tomorrow,
    # Adjust days to flower for three plants that were recorded two days late (but were flagged as 'flowered_yesterday').
    days_to_flower = ifelse( grepl("2 days ago,", comment), flowered_yesterday-1, days_to_flower)

  ) %>%
  filter(
    ! did_not_germinate, ! did_not_flower
  ) %>%
  select(
    id, replicate, tray, position, genotype, generation, cohort, days_to_flower
  )

flowering_time <- flowering_time %>%
  # Swap some names that I think were mixed up after validating genotypes.
  mutate(
    # genotype = case_when(
    #   genotype == "175x6913 rep1"  ~ "8307x9481 rep1",
    #   genotype == "175x6913 rep2"  ~ "8307x9481 rep2",
    #   genotype == "8307x9481 rep1" ~ "175x6913 rep1",
    #   genotype == "8307x9481 rep2" ~ "175x6913 rep2",
    #   genotype == "9470x6107 rep1" ~ "9470x7002 rep1",
    #   .default = genotype
    # ),
    genotype = gsub(" ", "_", genotype)
  )
