library(tidyverse)
library(gamlss)
setwd(here::here())

sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_liu2021.csv")
sub_luo2017 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_luo2017.csv")
sub_intmd <- rbind(sub_luo2017, sub_liu2021) %>% dplyr::filter(cell_MF > 0 & cell_MF < 1)
sub_intmd %>% 
  ggplot(aes(as.factor(med_cov), cell_MF)) +
  geom_violin(alpha = 0.5)
ggsave("plots/estim_emiBetaPrior_inflatedBM/violin_intmdMF_vs_MedCov.png", width = 10, height = 5)

sub_intmd %>% group_by(SubType, med_cov) %>% 
  summarise(med_cov = max(med_cov)) %>%
  View()

#################################################
##### fit beta mixture model on each subtype ####
#################################################

## install if not
# devtools::install_github("https://github.com/palmerimatthew/BetaMixture")

system.time(fit <- BetaMixture::BM_Fit(x = sub_intmd[SubType=="DG dg-all", cell_MF], 
                                       K = 2, threshold = 0.001, seed = 1))
fit

source("code/estim_emiBetaPrior_inflatedBM/helper_functions/BM_fit_modified.R")
system.time(fit <- BM_Fit(x = sub_intmd[SubType=="DG dg-all", cell_MF], 
                          threshold = 0.01, seed = 1))
fit


