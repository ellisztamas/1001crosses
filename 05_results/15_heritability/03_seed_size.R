#' Calculate broad-sense heritability for seed size
#'
#' This collates raw seed-size data, estimates variance explained by
#' genetics and experimental design, then expresses these estimates as a proportion
#' of total phenotypic variation on the scale of the data.
#' This is done separately for parents, and both cohorts of F9s.
#'
#' Tom Ellis, 16th September 2025


library('lme4')
library('QGglmm')

# Run formatting code used to create BLUPs.
# This saves a lot of boilerplate, but takes a long time because there is a
# heavy GLM to run, which is not actually not needed here.
source("03_processing/06_process_phenotypes/05_seed_size.R")

pops <- c("parent", "rep1", "rep2")

# Empty list to store variance components for the three groups
list_pve <- vector('list', 3)
names(list_pve) <- pops


for (group in pops){

  # Model flowering time as a function of genotype, cohort and tray
  mod_ss <- lmer(
    surface_area ~ (1 | rep/tray) + (1| genotype/label),
    data=raw_seed_sizes,
    subset = cohort==group
  )
  # Extract variance components on the expected scale
  vc <- as_tibble(VarCorr(mod_ss))
  vc$pve <- vc$vcov / sum(vc$vcov)

  # Extract variance components on the data scale
  # The model uses and identity link, so data and expected scales are identical
  vc <- as.data.frame(lme4::VarCorr(mod_ss))
  pve_data_scale <- vc$vcov / sum(vc$vcov)

  list_pve[[group]] <- pve_data_scale
}

# Data from with a row for the three groups, and columns for var explained by
# each component
h2 <- do.call('rbind', list_pve)
colnames(h2) <- c('Mother', 'Genotype', 'Tray', 'Replicate', 'Residuals')
h2 <- as_tibble(h2, rownames = "group")

# Plot a stacked bar chart with group on the x and var explained on the y
# Results are nearly identical for the three groups
# Genetics explains about 80% of the variance.
h2 %>%
  pivot_longer(Mother:Residuals, names_to = "component", values_to = "pve") %>%
  ggplot(aes(x = group, y = pve, fill=component)) +
  geom_col()+
  lims(
    y = c(0,1)
  ) +
  labs(
    y = "Prop. variance explained"
  ) +
  theme(
    axis.title.x = element_blank()
  )
