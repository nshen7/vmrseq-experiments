suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(data.table))
setwd(here::here())
set.seed(2022)

################################
#### Read in processed data ####
################################
sub_meth <- readRDS("data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_methClust.rds")

## !! Since the nu parameter in ZIBB for meth subpop is negligible, we just run BB regression

##############################################
##### BB regression for methylated subpop ####
##############################################

show.link("BB")
mod_m2 <- gamlss(y ~ log(med_cov), data = sub_meth, 
                 sigma.formula = ~ log(med_cov),
                 family=BB(mu.link = "logit", sigma.link = "log"),
                 control = gamlss.control(n.cyc = 200))
saveRDS(mod_m2, file = "code/estim_emiBetaPrior_ZIBBregression/model_meth_BBregression.rds")

summary(mod_m2)
# ******************************************************************
#   Family:  c("BB", "Beta Binomial") 
# 
# Call:  gamlss(formula = y ~ log(med_cov), sigma.formula = ~log(med_cov),  
#               family = BB(mu.link = "logit", sigma.link = "log"),  
#               data = sub_meth, control = gamlss.control(n.cyc = 200)) 
# 
# Fitting method: RS() 
# 
# ------------------------------------------------------------------
#   Mu link function:  logit
# Mu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.805278   0.007357  381.33   <2e-16 ***
#   log(med_cov) -0.088402   0.002331  -37.92   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Sigma link function:  log
# Sigma Coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)  -2.303557   0.015621 -147.470  < 2e-16 ***
#   log(med_cov) -0.026194   0.004774   -5.487  4.1e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   No. of observations in the fit:  1116502 
# Degrees of Freedom for the fit:  4
# Residual Deg. of Freedom:  1116498 
# at cycle:  10 
# 
# Global Deviance:     3634069 
# AIC:     3634077 
# SBC:     3634124 
# ******************************************************************
  
## diagnostics
png("plots/estim_emiBetaPrior_ZIBBregression/diagnosticPlots_methPop_BBregression.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_m2)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  -7.407699e-05 
# variance   =  0.9857013 
# coef. of skewness  =  -0.111258 
# coef. of kurtosis  =  3.151973 
# Filliben correlation coefficient  =  0.998386 
# ******************************************************************
  
sub_meth_smr <- sub_meth %>% 
  mutate(mu_fit = fitted(mod_m2, "mu"),
         sigma_fit = fitted(mod_m2, "sigma")) %>%
  group_by(DataSource, CellClass, SubType, med_cov) %>% 
  summarise(mean = mean(cell_MF), 
            var = var(cell_MF),
            mu_fit = max(mu_fit),
            sigma_fit = max(sigma_fit))

sub_meth_smr %>% ggplot(aes(mean, mu_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.91, 0.95) + ylim(0.91, 0.95) + xlab("MoM estimate of mu from IB") + ylab("1-Estimate of mu from BB regression")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_methPop_BBreg_vs_BetaMoM_paramMu.png", width = 5, height = 5)
sub_meth_smr %>% ggplot(aes(sqrt(var/mean/(1-mean)), sigma_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.089, 0.45) + ylim(0.089, 0.45) + xlab("MoM estimate of sigma from IB") + ylab("Estimate of sigma from BB regression")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_methPop_BBreg_vs_BetaMoM_paramSigma.png", width = 5, height = 5)

## regression estimates vs. median coverage
sub_meth_smr %>% ggplot() + 
  geom_point(aes(log(med_cov), mean), color = "darkgrey") + 
  geom_point(aes(log(med_cov), mu_fit), color = "blue") + 
  xlab("Log median across-cell coevrage") + ylab("Ests of mu from ZIBB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_methPop_paramMu_vs_medCov_BBreg&BetaMoM.png", width = 5, height = 5)
sub_meth_smr %>% ggplot() + 
  geom_point(aes(log(med_cov), sqrt(var/mean/(1-mean))), color = "darkgrey") + 
  geom_point(aes(log(med_cov), sigma_fit), color = "blue") + 
  xlab("Log median across-cell coevrage") + ylab("Estimates of sigma from ZIBB (blue) and MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_methPop_paramSigma_vs_medCov_BBreg&BetaMoM.png", width = 5, height = 5)



## Inspectiing null model
mod_m2_0 <- gamlss(y ~ 1, data = sub_meth, 
                   sigma.formula = ~ 1, 
                   family = BB,
                   control = gamlss.control(n.cyc = 200))
summary(mod_m2_0)
# Global Deviance:     3635965 
# AIC:     3635969 
# SBC:     3635993 

png("plots/estim_emiBetaPrior_ZIBBregression/diagnosticPlots_methPop_BBregression_nullModel.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_m2_0)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  0.002022558 
# variance   =  0.9815087 
# coef. of skewness  =  -0.1220704 
# coef. of kurtosis  =  3.184832 
# Filliben correlation coefficient  =  0.9982051 
# ******************************************************************
  