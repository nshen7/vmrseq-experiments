library(tidyverse)
library(gamlss)
library(gamlss.tr)
setwd(here::here())

## generate functions for truncated BEOI and BEZI distributions
gen.trun(par = 0.5, family = "BEZI", type = "right")
# A truncated family of distributions from BEZI has been generated 
# and saved under the names:  
#   dBEZItr pBEZItr qBEZItr rBEZItr BEZItr 
# The type of truncation is right 
# and the truncation parameter is 0.5  
gen.trun(par = 0.5, family = "BEOI", type = "left")
# A truncated family of distributions from BEOI has been generated 
# and saved under the names:  
#   dBEOItr pBEOItr qBEOItr rBEOItr BEOItr 
# The type of truncation is left 
# and the truncation parameter is 0.5 

########################################
# ==== fit inflated-beta regression ====
########################################
sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_liu2021.csv")
sub_luo2017 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_luo2017.csv")
sub_unmeth <- rbind(sub_luo2017, sub_liu2021) %>% dplyr::filter(cell_MF < 0.5)
sub_meth <- rbind(sub_luo2017, sub_liu2021) %>% dplyr::filter(cell_MF > 0.5)

sub_unmeth %>% filter(cell_MF > 0) %>% 
  ggplot(aes(as.factor(med_cov), cell_MF)) +
  geom_jitter(width = 0.2, height = 0.01, alpha = 0.3, size = 0.5, color = "steelblue") +
  geom_violin(alpha = 0.5)
ggsave("plots/estim_emiBetaPrior_IBregression/violin_unmethPop_non0MF_vs_MedCov.png", width = 7, height = 5)

sub_meth %>% filter(cell_MF < 1) %>% 
  ggplot(aes(as.factor(med_cov), cell_MF)) +
  geom_jitter(width = 0.2, height = 0.01, alpha = 0.1, size = 0.5, color = "steelblue") + 
  geom_violin(alpha = 0.5)
ggsave("plots/estim_emiBetaPrior_IBregression/violin_methPop_non1MF_vs_MedCov.png", width = 7, height = 5)

#######################################
##### regression for unmethylated prior

show.link(family = "BEZItr")
mod_u <- gamlss(cell_MF ~ log(med_cov), data = sub_unmeth[sample(1:nrow(sub_unmeth), 20000)], 
                sigma.formula = ~ log(med_cov), nu.formula = ~ log(med_cov), 
                family=BEZItr(nu.link = "identity", mu.link = "identity", sigma.link = "inverse"), 
                control = gamlss.control(n.cyc = 40, c.crit = 0.01))

## diagnostic plots
summary(mod_u)
png("plots/estim_emiBetaPrior_IBregression/diagnosticPlots_unmethPop_IBregressgion.png", width = 1000, height = 700)
plot(mod_u)
dev.off()

sub_unmeth_smr <- sub_unmeth %>% 
  mutate(nu_fit = fitted(mod_u, "nu"),
         mu_fit = fitted(mod_u, "mu"),
         sigma_fit = fitted(mod_u, "sigma")) %>%
  group_by(DataSource, CellClass, SubType, med_cov) %>% 
  summarise(weight = sum(cell_MF==0) / length(cell_MF), 
            mean = mean(cell_MF[cell_MF>0]), 
            var = var(cell_MF[cell_MF>0]),
            nu_fit = max(nu_fit),
            mu_fit = max(mu_fit),
            sigma_fit = max(sigma_fit))
# MOM estimates (plug in mean and var) vs. regression estimates
sub_unmeth_smr %>% ggplot(aes(weight, nu_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.3, 0.7) + ylim(0.3, 0.7) + xlab("MOM estimate of nu") + ylab("Estimate of nu from inflated beta regression")
ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_IBreg_vs_MOM_paramNu.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot(aes(mean, mu_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.1, 0.3) + ylim(0.1, 0.3) + xlab("MOM estimate of mu") + ylab("Estimate of mu from inflated beta regression")
ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_IBreg_vs_MOM_paramMu.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot(aes(mu_fit*(1-mu_fit)/var - 1, sigma_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(4, 13) + ylim(4, 13) + xlab("MOM estimate of sigma") + ylab("Estimate of sigma from inflated beta regression")
ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_IBreg_vs_MOM_paramSigma.png", width = 5, height = 5)

#####################################
##### regression for methylated prior

show.link(family = "BEOI")
mod_m <- gamlss(cell_MF ~ log(med_cov), data = sub_meth, 
                sigma.formula = ~ log(med_cov), nu.formula = ~ log(med_cov), 
                family=BEOI(nu.link = "identity", mu.link = "identity", sigma.link = "inverse"))
summary(mod_m)

png("plots/estim_emiBetaPrior_IBregression/diagnosticPlots_methPop_IBregressgion.png", width = 1000, height = 700)
plot(mod_m)
dev.off()

sub_meth_smr <- sub_meth %>% 
  mutate(nu_fit = fitted(mod_m, "nu"),
         mu_fit = fitted(mod_m, "mu"),
         sigma_fit = fitted(mod_m, "sigma")) %>%
  group_by(DataSource, CellClass, SubType, med_cov) %>% 
  summarise(weight = sum(cell_MF==1) / length(cell_MF), 
            mean = mean(cell_MF[cell_MF<1]), 
            var = var(cell_MF[cell_MF<1]),
            nu_fit = max(nu_fit),
            mu_fit = max(mu_fit),
            sigma_fit = max(sigma_fit))
# MOM estimates (plug in mean and var) vs. regression estimates
sub_meth_smr %>% ggplot(aes(weight, nu_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlim(0.1, 0.7) + ylim(0.1, 0.7) + xlab("MOM estimate of nu") + ylab("Estimate of nu from inflated beta regression")
ggsave("plots/estim_emiBetaPrior_IBregression/point_methPop_IBreg_vs_MOM_paramNu.png", width = 5, height = 5)
sub_meth_smr %>% ggplot(aes(mean, mu_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlim(0.8, 0.95) + ylim(0.8, 0.95) + xlab("MOM estimate of nu") + ylab("Estimate of nu from inflated beta regression")
ggsave("plots/estim_emiBetaPrior_IBregression/point_methPop_IBreg_vs_MOM_paramMu.png", width = 5, height = 5)
sub_meth_smr %>% ggplot(aes(mu_fit*(1-mu_fit)/var - 1, sigma_fit, color = med_cov)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0, color = "red") +
  xlim(7, 16) + ylim(7, 16) + xlab("MOM estimate of nu") + ylab("Estimate of nu from inflated beta regression")
ggsave("plots/estim_emiBetaPrior_IBregression/point_methPop_IBreg_vs_MOM_paramSigma.png", width = 5, height = 5)

