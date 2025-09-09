#' Calculate broad-sense heritability for flowering time
#'
#' This collates raw flowering time data, estimates variance explained by
#' genetics and experimental design, then expresses these estimates as a proportion
#' of total phenotypic variation on the scale of the data data.


library('lme4')
library('QGglmm')

# Collate and format raw flowering-time data
source("03_processing/06_process_phenotypes/02_rosette_size.R")

pops <- c("parent", "rep1", "rep2")

# Empty list to store variance components for the three groups
list_pve <- vector('list', 3)
names(list_pve) <- pops

# Variance components for three cohorts
for (group in pops){

  # Model flowering time as a function of genotype, cohort and tray
  rosette_model <- lmer(
    sqrt(rosette) ~ (1 | tray) + (1 | genotype) + (1|id),
    data = raw_rosette_size, subset = cohort==group
  )
  # Extract variance components on the data scale
  # The model uses and identity link, so data and expected scales are identical
  vc <- as.data.frame(lme4::VarCorr(rosette_model))
  pve_data_scale <- vc$vcov / sum(vc$vcov)

  list_pve[[group]] <- pve_data_scale
}

# Data from with a row for the three groups, and columns for var explained by
# each component
h2 <- do.call('rbind', list_pve)
colnames(h2) <- c('Measurement', 'Genotype', 'Tray', 'residuals')
h2 <- as_tibble(h2, rownames = "group")

# Plot a stacked bar chart with group on the x and var explained on the y
# Results are nearly identical for the three groups
# Genetics explains about 80% of the variance.
h2 %>%
  pivot_longer(Measurement:Tray, names_to = "component", values_to = "pve") %>%
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
