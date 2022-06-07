library(tidyverse)
library(gamlss)
setwd(here::here())

sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_liu2021.csv")
sub_luo2017 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_luo2017.csv")
sub_unmeth <- rbind(sub_luo2017, sub_liu2021) %>% dplyr::filter(cell_MF < 0.5 & cell_MF > 0)
sub_meth <- rbind(sub_luo2017, sub_liu2021) %>% dplyr::filter(cell_MF > 0.5 & cell_MF < 1)

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
show.link(family = "BE")
mod_u <- gamlss(cell_MF ~ log(med_cov), data = sub_unmeth, 
                sigma.formula = ~ log(med_cov), 
                family=BE(mu.link = "identity", sigma.link = "inverse"))
## diagnostic plots
summary(mod_u)
png("plots/estim_emiBetaPrior_IBregression/diagnosticPlots_unmethPop_BEregression.png", width = 1000, height = 700)
plot(mod_u)
dev.off()

wp(mod_u, ylim.all = 1)

sub_unmeth_smr <- sub_unmeth %>% 
  mutate(mu_fit = fitted(mod_u, "mu"),
         sigma_fit = fitted(mod_u, "sigma")) %>%
  group_by(DataSource, CellClass, SubType, med_cov) %>% 
  summarise(mean = mean(cell_MF[cell_MF>0]), 
            var = var(cell_MF[cell_MF>0]),
            mu_fit = max(mu_fit),
            sigma_fit = max(sigma_fit),
            alpha_fit = mu_fit*(1/sigma_fit^2-1),
            beta_fit = (1-mu_fit)*(1/sigma_fit^2-1))
# MOM estimates (plug in mean and var) vs. regression estimates
sub_unmeth_smr %>% ggplot(aes(mean, mu_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.1, 0.3) + ylim(0.1, 0.3) + xlab("MOM estimate of mu") + ylab("Estimate of mu from inflated beta regression")
# ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_IBreg_vs_MOM_paramMu.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot(aes(sqrt(var/mean/(1-mean)), sigma_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.25, 0.5) + ylim(0.25, 0.5) + xlab("MOM estimate of sigma") + ylab("Estimate of sigma from inflated beta regression")
# ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_IBreg_vs_MOM_paramSigma.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot(aes((mean*(1-mean)/var - 1) * mean, alpha_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.5, 3.2) + ylim(0.5, 3.2) + xlab("MOM estimate of alpha") + ylab("Estimate of alpha from inflated beta regression")
sub_unmeth_smr %>% ggplot(aes((mean*(1-mean)/var - 1) * (1-mean), beta_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(3,9) + ylim(3,9) + xlab("MOM estimate of alpha") + ylab("Estimate of alpha from inflated beta regression")

#####################################
##### regression for methylated prior
sub_meth_smr0 <- sub_meth %>% 
  group_by(DataSource, CellClass, SubType, med_cov) %>% 
  summarise(weight = sum(cell_MF==1) / length(cell_MF), 
            mean = mean(cell_MF), 
            var = var(cell_MF))

sub_meth_smr0 %>% 
  ggplot(aes(log(med_cov), sqrt(var/mean/(1-mean))^(-1))) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  ylab("Inverse MOM estimate of sigma") + xlab("Log(median coverage)")

show.link(family = "BE")
mod_m <- gamlss(cell_MF ~ log(med_cov), data = sub_meth, 
                sigma.formula = ~ log(med_cov), 
                family=BE(mu.link = "identity", sigma.link = "inverse"))
summary(mod_m)

png("plots/estim_emiBetaPrior_IBregression/diagnosticPlots_methPop_BEregression.png", width = 1000, height = 700)
plot(mod_m)
dev.off()

term.plot(mod_m, "mu")
term.plot(mod_m, "sigma")
wp(mod_m, ylim.all = 1)

sub_meth_smr <- sub_meth %>% 
  mutate(mu_fit = fitted(mod_m, "mu"),
         sigma_fit = fitted(mod_m, "sigma")) %>%
  group_by(DataSource, CellClass, SubType, med_cov) %>% 
  summarise(mean = mean(cell_MF), 
            var = var(cell_MF),
            mu_fit = max(mu_fit),
            sigma_fit = max(sigma_fit),
            alpha_fit = mu_fit*(1/sigma_fit^2-1),
            beta_fit = (1-mu_fit)*(1/sigma_fit^2-1))

# MOM estimates (plug in mean and var) vs. regression estimates
sub_meth_smr %>% ggplot(aes(mean, mu_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlim(0.8, 0.95) + ylim(0.8, 0.95) + xlab("MOM estimate of mu") + ylab("Estimate of mu from inflated beta regression")
sub_meth_smr %>% ggplot(aes(sqrt(var/mean/(1-mean)), sigma_fit, color = med_cov)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0, color = "red") +
  xlim(0.23, 0.35) + ylim(0.23, 0.35) + xlab("MOM estimate of sigma") + ylab("Estimate of sigma from inflated beta regression")
sub_meth_smr %>% ggplot(aes((mean*(1-mean)/var - 1) * mean, alpha_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(7, 13) + ylim(7, 13) + xlab("MOM estimate of alpha") + ylab("Estimate of alpha from inflated beta regression")
sub_meth_smr %>% ggplot(aes((mean*(1-mean)/var - 1) * (1-mean), beta_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.7, 3.0) + ylim(0.7, 3.0) + xlab("MOM estimate of alpha") + ylab("Estimate of alpha from inflated beta regression")

sub_meth_smr %>% ggplot(aes(mu_fit, sigma_fit, color = med_cov)) + 
  geom_point() 