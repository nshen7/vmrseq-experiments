suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(data.table))
setwd(here::here())
set.seed(2022)


##################################################
##### ZIBB regression for methylated subpop ####
##################################################

## Read in training data
sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_liu2021_clusterMeth.csv")
sub_luo2017m <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017mice_clusterMeth.csv")
sub_luo2017h <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017human_clusterMeth.csv")
sub_meth <- rbind(sub_liu2021, sub_luo2017m, sub_luo2017h) %>% 
  as_tibble() %>%
  mutate(y = as.matrix(data.frame(cell_cov - cell_meth, cell_meth))) %>% 
  # filter(log(med_cov) <= 4) %>% 
  # slice_sample(n = 20000) %>% 
  dplyr::select(-CellClass) # CellClass contains NA from liu2017 human data so removed

## run ZIBB regression
show.link("ZIBB")
# mod_m <- gamlss(y ~ lo(~log(med_cov), degree = 1), data = sub_meth,
#                 sigma.formula = ~ 1,
#                 nu.formula = ~ lo(~log(med_cov), degree = 1),
#                 family = ZIBB(nu.link = "logit", mu.link = "logit", sigma.link = "log"),
#                 n.cyc = 50, method = mixed(), c.crit = 0.01) # n = 50000: 3832825.pbsha.*; all: 3840383.pbsha.ib.sockeye

# mod_m <- gamlss(y ~ pb(log(med_cov), inter = 3, degree = 2), data = sub_meth,
#                 sigma.formula = ~ pb(log(med_cov), inter = 3, degree = 2),
#                 nu.formula = ~ pb(log(med_cov), inter = 3, degree = 2),
#                 family = ZIBB(nu.link = "logit", mu.link = "logit", sigma.link = "log"),
#                 n.cyc = 100)
mod_m <- gamlss(y ~ 1, data = sub_meth,
                sigma.formula = ~ 1,
                nu.formula = ~ 1,
                family = ZIBB(nu.link = "logit", mu.link = "logit", sigma.link = "log"),
                n.cyc = 100)
# 2 hrs: 3866978.pbsha.ib.sockeye
# 120 hrs: 3866979.pbsha.ib.sockeye
saveRDS(mod_m, file = "code/estim_emiBetaPrior_ZIBBregression2/model_meth_ZIBBregression.rds")

## diagnostics
summary(mod_m)


# png("plots/estim_emiBetaPrior_ZIBBregression2/diagnosticPlots_methPop_ZIBBregression.png", width = 1000, height = 700)
# par(mar=c(2,1.5,1.5,1.5)); plot(mod_m)
# dev.off()


sub_meth_smr <- sub_meth %>%
  mutate(nu_fit = fitted(mod_m, "nu"),
         mu_fit = fitted(mod_m, "mu"),
         sigma_fit = fitted(mod_m, "sigma")) %>%
  group_by(DataSource, SubType, med_cov) %>%
  summarise(nu_mle = max(mle_nu),
            mu_mle = max(mle_mu),
            sigma_mle = max(mle_sigma),
            nu_fit = max(nu_fit),
            mu_fit = max(mu_fit),
            sigma_fit = max(sigma_fit))

## regression estimates vs. median coverage
sub_meth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), nu_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), nu_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Ests of nu from ZIBB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_methPop_paramNu_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)
sub_meth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), mu_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), mu_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Ests of mu from ZIBB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_methPop_paramMu_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)
sub_meth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), sigma_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), sigma_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Estimates of sigma from ZIBB (blue) and MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_methPop_paramSigma_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)

## Inspectiing null model
# mod_m0 <- gamlss(y ~ 1, data = sub_meth, 
#                  sigma.formula = ~ 1, nu.formula = ~ 1, 
#                  family=ZIBB,
#                  control = gamlss.control(n.cyc = 200))
# summary(mod_m0)
# 
# png("plots/estim_emiBetaPrior_ZIBBregression2/diagnosticPlots_methPop_ZIBBregression_nullModel.png", width = 1000, height = 700)
# par(mar=c(2,1.5,1.5,1.5)); plot(mod_m0)
# dev.off()
