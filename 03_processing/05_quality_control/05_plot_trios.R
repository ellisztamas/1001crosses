library("tidyverse")

trio_files <- Sys.glob("03_processing/05_quality_control/output/data_to_check_trios/*_check_trio.csv")

i <- 1
trio <- read_csv(trio_files[i])
trio %>%
    pivot_longer(progeny_vs_p1:p1_vs_p2) %>%
    # filter(name != "progeny_vs_p2") %>%
    ggplot(aes(x = start, y = value, colour = name)) +
    geom_line() +
    # labs(
    #   title = basename(trio_files[i])
    # ) +
    facet_grid(rows=vars(chr))

hist(trio$progeny_vs_p1[1:10])

rt_files <- Sys.glob("03_processing/05_quality_control/output/random_trios/*csv")
rt_files[1] %>%
  read_csv() %>%
  pivot_longer(progeny_vs_p1:p1_vs_p2) %>%
  ggplot(aes(x = start, y = value, colour = name)) +
  geom_line() +
  facet_grid(rows=vars(chr))



rank_files <- Sys.glob(
  "03_processing/05_quality_control/output/rank_candidates/*csv"
  )

i <- 4
(rank <- read_csv(rank_files[1]))
exp(rank$score)
which(rank$putative_parents)
which(rank$parent2 == 7424)
7424 %in% rank$parent2
