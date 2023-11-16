suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(data.table))
setwd(here::here())
set.seed(2022)

###################################################
##### ZIBB regression for unmethylated subpop #####
###################################################

## Read in training data
sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_liu2021_clusterUnmeth.csv")
sub_luo2017m <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017mice_clusterUnmeth.csv")
sub_luo2017h <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017human_clusterUnmeth.csv")
sub_unmeth <- rbind(sub_liu2021, sub_luo2017m, sub_luo2017h) %>% 
  as_tibble() %>%
  mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth))) %>% 
  # slice_sample(n = 20000) %>% # filter(log(med_cov) <= 4.5) %>%
  dplyr::select(-CellClass) # CellClass contains NA from liu2017 human data so removed

## Fit ZIBB regression on nu, mu, sigma against log(med_cov)
show.link("ZIBB")
mod_u <- gamlss(y ~ cs(log(med_cov)), data = sub_unmeth,
                sigma.formula = ~ cs(log(med_cov)),
                nu.formula = ~ cs(log(med_cov)),
                family = ZIBB(nu.link = "logit", mu.link = "logit", sigma.link = "log"),
                n.cyc = 100)
# 2 hrs: 3866974.pbsha.ib.sockeye
# 120 hrs: 3866976.pbsha.ib.sockeye
saveRDS(mod_u, file = "code/estim_emiBetaPrior_ZIBBregression2/model_unmeth_ZIBBregression.rds")

## diagnostics
summary(mod_u)


# png("plots/estim_emiBetaPrior_ZIBBregression2/diagnosticPlots_unmeth_ZIBBregression.png", width = 1000, height = 700)
# par(mar=c(2,1.5,1.5,1.5)); plot(mod_u)
# dev.off()

sub_unmeth_smr <- sub_unmeth %>%
  mutate(nu_fit = fitted(mod_u, "nu"),
         mu_fit = fitted(mod_u, "mu"),
         sigma_fit = fitted(mod_u, "sigma")) %>%
  group_by(DataSource, SubType, med_cov) %>%
  summarise(nu_mle = max(mle_nu),
            mu_mle = max(mle_mu),
            sigma_mle = max(mle_sigma),
            nu_fit = max(nu_fit),
            mu_fit = max(mu_fit),
            sigma_fit = max(sigma_fit))

## regression estimates vs. median coverage
sub_unmeth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), nu_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), nu_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Ests of nu from ZIBB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_unmeth_paramNu_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), mu_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), mu_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Ests of mu from ZIBB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_unmeth_paramMu_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), sigma_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), sigma_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Estimates of sigma from ZIBB (blue) and MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_unmeth_paramSigma_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)


## Inspecting null model
mod_u0 <- gamlss(y ~ 1, data = sub_unmeth,
                 sigma.formula = ~ 1, nu.formula = ~ 1,
                 family=ZIBB,
                 control = gamlss.control(n.cyc = 100))
summary(mod_u0)

png("plots/estim_emiBetaPrior_ZIBBregression2/diagnosticPlots_unmeth_ZIBBregression_nullModel.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_u0)
dev.off()


