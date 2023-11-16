suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(data.table))
setwd(here::here())
set.seed(2022)

################################
#### Read in processed data ####
################################
sub_unmeth <- readRDS("data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_unmethClust.rds")

##################################################
##### ZIBB regression for unmethylated subpop ####
##################################################

show.link("ZIBB")
mod_u <- gamlss(y ~ log(med_cov), data = sub_unmeth,
                sigma.formula = ~ log(med_cov), nu.formula = ~ log(med_cov),
                family=ZIBB(nu.link = "logit", mu.link = "logit", sigma.link = "log"),
                control = gamlss.control(n.cyc = 200))
saveRDS(mod_u, file = "code/estim_emiBetaPrior_ZIBBregression/model_unmeth_ZIBBregression.rds")
## diagnostics
summary(mod_u)
# ******************************************************************
#   Family:  c("ZIBB", "Zero Inflated Beta Binomial") 
# 
# Call:  gamlss(formula = y ~ log(med_cov), sigma.formula = ~log(med_cov),  
#               nu.formula = ~log(med_cov), family = ZIBB(nu.link = "logit",  
#                                                         mu.link = "logit", sigma.link = "log"), data = sub_unmeth,  
#               control = gamlss.control(n.cyc = 200)) 
# 
# Fitting method: RS() 
# 
# ------------------------------------------------------------------
#   Mu link function:  logit
# Mu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.696779   0.023408  -29.77   <2e-16 ***
#   log(med_cov) -0.348365   0.007497  -46.47   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Sigma link function:  log
# Sigma Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -3.53347    0.04504  -78.45   <2e-16 ***
#   log(med_cov)  0.59335    0.01305   45.46   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Nu link function:  logit 
# Nu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.76683    0.04287   41.22   <2e-16 ***
#   log(med_cov) -0.82164    0.01574  -52.21   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   No. of observations in the fit:  226120 
# Degrees of Freedom for the fit:  6
# Residual Deg. of Freedom:  226114 
# at cycle:  15 
# 
# Global Deviance:     804496.1 
# AIC:     804508.1 
# SBC:     804570.1 
# ******************************************************************
png("plots/estim_emiBetaPrior_ZIBBregression/diagnosticPlots_unmethPop_ZIBBregression.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_u)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  0.0051781 
# variance   =  0.9836178 
# coef. of skewness  =  -0.1109042 
# coef. of kurtosis  =  2.572494 
# Filliben correlation coefficient  =  0.9937089 
# ******************************************************************
  
ggplot() +
  geom_violin(aes(as.factor(sub_unmeth$med_cov), mod_u$residuals))

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
# MoM estimates (plug in mean and var) vs. regression estimates
sub_unmeth_smr %>% ggplot(aes(weight, nu_fit, color = med_cov)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlim(0.15, 0.66) + ylim(0.15, 0.66) + xlab("MoM estimate of nu from IB") + ylab("Estimate of nu from ZIBB regression")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_unmethPop_ZIBBreg_vs_IBMoM_paramNu.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot(aes(mean, mu_fit, color = med_cov)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlim(0.1, 0.26) + ylim(0.1, 0.26) + xlab("MoM estimate of mu from IB") + ylab("Estimate of mu from ZIBB regression")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_unmethPop_ZIBBreg_vs_IBMoM_paramMu.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot(aes(sqrt(var/mean/(1-mean)), sigma_fit, color = med_cov)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlim(0.1, 0.45) + ylim(0.1, 0.45) + xlab("MoM estimate of sigma from IB") + ylab("Estimate of sigma from ZIBB regression")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_unmethPop_ZIBBreg_vs_IBMoM_paramSigma.png", width = 5, height = 5)

## regression estimates vs. median coverage
sub_unmeth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), weight), color = "darkgrey") +
  geom_point(aes(log(med_cov), nu_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Ests of nu from ZIBB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_unmethPop_paramNu_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), mean), color = "darkgrey") +
  geom_point(aes(log(med_cov), mu_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Ests of mu from ZIBB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_unmethPop_paramMu_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), sqrt(var/mean/(1-mean))), color = "darkgrey") +
  geom_point(aes(log(med_cov), sigma_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Estimates of sigma from ZIBB (blue) and MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_unmethPop_paramSigma_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)


## Inspecting null model
mod_u0 <- gamlss(y ~ 1, data = sub_unmeth,
                 sigma.formula = ~ 1, nu.formula = ~ 1,
                 family=ZIBB,
                 control = gamlss.control(n.cyc = 200))
summary(mod_u0)
# Global Deviance:     807493.4 
# AIC:     807499.4 
# SBC:     807530.3 

png("plots/estim_emiBetaPrior_ZIBBregression/diagnosticPlots_unmethPop_ZIBBregression_nullModel.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_u0)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  0.00541673 
# variance   =  0.9960295 
# coef. of skewness  =  -0.1332435 
# coef. of kurtosis  =  2.545661 
# Filliben correlation coefficient  =  0.992856 
# ******************************************************************
  

