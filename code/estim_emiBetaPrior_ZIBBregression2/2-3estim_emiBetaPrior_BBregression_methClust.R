suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(data.table))
setwd(here::here())
set.seed(2022)

## !! Since the nu parameter in ZIBB for meth subpop is negligible, we just run BB regression

##############################################
##### BB regression for methylated subpop ####
##############################################

## Read in training data
sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_liu2021_clusterMeth.csv")
sub_luo2017m <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017mice_clusterMeth.csv")
sub_luo2017h <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017human_clusterMeth.csv")
sub_meth <- rbind(sub_liu2021, sub_luo2017m, sub_luo2017h) %>% 
  as_tibble() %>%
  mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth))) %>% 
  # slice_sample(n = 20000) %>%
  # filter(log(med_cov) <= 4) %>%
  dplyr::select(-CellClass) # CellClass contains NA from liu2017 human data so removed

saveRDS(sub_meth, "data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_methClust.rds")

show.link("BB")
mod_m2 <- gamlss(y ~ 1, data = sub_meth,
                 sigma.formula = ~ 1, 
                 family = BB(mu.link = "logit", sigma.link = "log"),
                 n.cyc = 50) 
# 2 hrs: 3866982.pbsha.ib.sockeye
# 120 hrs: 3866983.pbsha.ib.sockeye
saveRDS(mod_m2, file = "code/estim_emiBetaPrior_ZIBBregression2/model_meth_BBregression.rds")

summary(mod_m2)

## diagnostics
# png("plots/estim_emiBetaPrior_ZIBBregression2/diagnosticPlots_methPop_BBregression.png", width = 1000, height = 700)
# par(mar=c(2,1.5,1.5,1.5)); plot(mod_m2)
# dev.off()

sub_meth_smr <- sub_meth %>%
  mutate(mu_fit = fitted(mod_m2, "mu"),
         sigma_fit = fitted(mod_m2, "sigma")) %>%
  group_by(DataSource, SubType, med_cov) %>%
  summarise(mu_mle = max(bbmle_mu),
            sigma_mle = max(bbmle_sigma),
            mu_fit = max(mu_fit),
            sigma_fit = max(sigma_fit))


## regression estimates vs. median coverage
sub_meth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), mu_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), mu_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Ests of mu from BB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_methPop_paramMu_vs_medCov_BBreg&IBMoM.png", width = 5, height = 5)
sub_meth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), sigma_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), sigma_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Estimates of sigma from BB (blue) and MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_methPop_paramSigma_vs_medCov_BBreg&IBMoM.png", width = 5, height = 5)



## Inspectiing null model
mod_m2_0 <- gamlss(y ~ 1, data = sub_meth, 
                   sigma.formula = ~ 1, 
                   family = BB,
                   control = gamlss.control(n.cyc = 200))
summary(mod_m2_0)


png("plots/estim_emiBetaPrior_ZIBBregression2/diagnosticPlots_methPop_BBregression_nullModel.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_m2_0)
dev.off()
