library(tidyverse)
library(gamlss)
library(gamlss.tr)
setwd(here::here())

sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_liu2021.csv")
sub_luo2017 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_luo2017.csv")
set.seed(2022)
sub_unmeth <- rbind(sub_luo2017, sub_liu2021) %>% 
  dplyr::filter(cell_MF < 0.5 & cell_MF > 0)
sub_meth <- rbind(sub_luo2017, sub_liu2021) %>% 
  dplyr::filter(cell_MF > 0.5 & cell_MF < 1)

#######################################
##### regression for unmethylated prior
gen.trun(par = 0.5, family = "BE", type = "right")
mod_u <- gamlss(cell_MF ~ log(med_cov), data = sub_unmeth, 
                sigma.formula = ~ log(med_cov), 
                family=BEtr(mu.link = "identity", sigma.link = "inverse"),
                control = gamlss.control(n.cyc = 200, mu.step = 1, sigma.step = 1))
## diagnostics
summary(mod_u)
png("plots/estim_emiBetaPrior_IBregression/diagnosticPlots_unmethPop_BEregression_truncated.png", width = 1000, height = 700)
plot(mod_u)
dev.off()
# ******************************************************************
#   Summary of the Quantile Residuals
# mean   =  -0.009001585 
# variance   =  0.9905133 
# coef. of skewness  =  0.2929081 
# coef. of kurtosis  =  1.991388 
# Filliben correlation coefficient  =  0.979232 
# ******************************************************************
hist(mod_u$residuals[sub_unmeth$med_cov==8], breaks = 20)
hist(mod_u$residuals[sub_unmeth$med_cov==66], breaks = 20)

ggplot() + 
  geom_violin(aes(as.factor(sub_unmeth$med_cov), resid(mod_u))) +
  xlab("Median coverage") + ylab("Quantile residuals")
ggsave("plots/estim_emiBetaPrior_IBregression/violin_unmethPop_BEtruncated_QRes_vs_medCov.png", width = 5, height = 5)

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
ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_BEtruncated_vs_MOM_paramMu.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot(aes(sqrt(var/mean/(1-mean)), sigma_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.25, 0.56) + ylim(0.25, 0.56) + xlab("MOM estimate of sigma") + ylab("Estimate of sigma from inflated beta regression")
ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_BEtruncated_vs_MOM_paramSigma.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot(aes((mean*(1-mean)/var - 1) * mean, alpha_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.5, 3.2) + ylim(0.5, 3.2) + xlab("MOM estimate of alpha") + ylab("Estimate of alpha from inflated beta regression")
sub_unmeth_smr %>% ggplot(aes((mean*(1-mean)/var - 1) * (1-mean), beta_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(3,9) + ylim(3,9) + xlab("MOM estimate of alpha") + ylab("Estimate of alpha from inflated beta regression")

#####################################
##### regression for methylated prior

gen.trun(par = 0.5, family = "BE", type = "left")
mod_m <- gamlss(cell_MF ~ log(med_cov), data = sub_meth, 
                sigma.formula = ~ log(med_cov), 
                family=BEtr(mu.link = "identity", sigma.link = "inverse"),
                control = gamlss.control(n.cyc = 40, mu.step = 1, sigma.step = 1))
summary(mod_m)

png("plots/estim_emiBetaPrior_IBregression/diagnosticPlots_methPop_BEregression_truncated.png", width = 1000, height = 700)
plot(mod_m)
dev.off()
# ******************************************************************
#   Summary of the Quantile Residuals
# mean   =  0.01562877 
# variance   =  1.005504 
# coef. of skewness  =  -0.9574977 
# coef. of kurtosis  =  3.337744 
# Filliben correlation coefficient  =  0.9604264 
# ******************************************************************
  
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
ggsave("plots/estim_emiBetaPrior_IBregression/point_methPop_BEtruncated_vs_MOM_paramMu.png", width = 5, height = 5)
sub_meth_smr %>% ggplot(aes(sqrt(var/mean/(1-mean)), sigma_fit, color = med_cov)) + 
  geom_point() + geom_abline(slope = 1, intercept = 0, color = "red") +
  xlim(0.23, 0.35) + ylim(0.23, 0.35) + xlab("MOM estimate of sigma") + ylab("Estimate of sigma from inflated beta regression")
ggsave("plots/estim_emiBetaPrior_IBregression/point_methPop_BEtruncated_vs_MOM_paramSigma.png", width = 5, height = 5)
sub_meth_smr %>% ggplot(aes((mean*(1-mean)/var - 1) * mean, alpha_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(7, 13) + ylim(7, 13) + xlab("MOM estimate of alpha") + ylab("Estimate of alpha from inflated beta regression")
sub_meth_smr %>% ggplot(aes((mean*(1-mean)/var - 1) * (1-mean), beta_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.7, 3.0) + ylim(0.7, 3.0) + xlab("MOM estimate of alpha") + ylab("Estimate of alpha from inflated beta regression")

sub_meth_smr %>% ggplot(aes(mu_fit, sigma_fit, color = med_cov)) + 
  geom_point() 