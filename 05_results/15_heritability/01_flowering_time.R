#' Calculate broad-sense heritability for flowering time
#'
#' This collates raw flowering time data, estimates variance explained by
#' genetics and experimental design, then expresses these estimates as a proportion
#' of total phenotypic variation on the scale of the data.


library('lme4')
library('QGglmm')

# Collate and format raw flowering-time data
source("03_processing/06_process_phenotypes/03_flowering_time.R")

pops <- c("parents", "rep1", "rep2")

# Empty list to store variance components for the three groups
list_pve <- vector('list', 3)
names(list_pve) <- pops


for (group in pops){

  # Model flowering time as a function of genotype, cohort and tray
  mod_ft <- glmer(
    days_to_flower ~  (1 | genotype) + (1| replicate/tray),
    data=flowering_time, family = 'poisson',
    subset = cohort==group
  )
  # Extract variance components on the expected scale
  vc <- as_tibble(VarCorr(mod_ft))
  vc$pve <- vc$vcov / sum(vc$vcov)

  # Use QgGLMM to get variance components on the data scale (including Poisson noise)
  # Each QGparams output gives the variance explained by each component divided
  # by total phenotypic variance on the data scale.
  vc$pve_data_scale <- NA
  for( i in 1:nrow(vc) ){
    vc$pve_data_scale[i] <- QGparams(
      mu = fixef(mod_ft),
      var.a = vc$vcov[i],
      var.p = sum(vc$vcov),
      model = 'Poisson.log',
      verbose = FALSE
    )$h2.obs
  }
  list_pve[[group]] <- vc$pve_data_scale
}

# Data from with a row for the three groups, and columns for var explained by
# each component
h2 <- do.call('rbind', list_pve )
colnames(h2) <- c('genotype', 'tray', 'replicate')
h2 <- as_tibble(h2, rownames = "group")

# Plot a stacked bar chart with group on the x and var explained on the y
# Results are nearly identical for the three groups
# Genetics explains about 80% of the variance.
h2 %>%
  pivot_longer(genotype:replicate, names_to = "component", values_to = "pve") %>%
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



