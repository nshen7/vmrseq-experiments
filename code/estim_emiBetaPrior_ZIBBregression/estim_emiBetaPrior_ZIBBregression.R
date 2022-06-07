library(tidyverse)
library(gamlss)
library(gamlss.tr)
setwd(here::here())

sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_liu2021.csv")
sub_luo2017 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_luo2017.csv")
sub_unmeth <- rbind(sub_luo2017, sub_liu2021) %>% 
  as_tibble() %>%
  dplyr::filter(cell_MF < 0.5) %>%
  mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth)))
sub_meth <- rbind(sub_luo2017, sub_liu2021) %>% 
  as_tibble() %>%
  dplyr::filter(cell_MF > 0.5) %>%
  mutate(y = as.matrix(data.frame(cell_cov - cell_meth, cell_meth)))

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
# (Intercept)  -0.77082    0.05272  -14.62   <2e-16 ***
#   log(med_cov) -0.32704    0.01680  -19.46   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Sigma link function:  log
# Sigma Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -3.47538    0.09981  -34.82   <2e-16 ***
#   log(med_cov)  0.58143    0.02891   20.11   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Nu link function:  logit 
# Nu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.67383    0.09791   17.10   <2e-16 ***
#   log(med_cov) -0.80262    0.03584  -22.39   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   No. of observations in the fit:  45370 
# Degrees of Freedom for the fit:  6
# Residual Deg. of Freedom:  45364 
# at cycle:  15 
# 
# Global Deviance:     161921.2 
# AIC:     161933.2 
# SBC:     161985.5 
# ******************************************************************

png("plots/estim_emiBetaPrior_ZIBBregression/diagnosticPlots_unmethPop_ZIBBregression.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_u)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  9.044864e-05 
# variance   =  0.9990006 
# coef. of skewness  =  -0.1242583 
# coef. of kurtosis  =  2.575818 
# Filliben correlation coefficient  =  0.9936238 
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
# Global Deviance:     163196.3 
# AIC:     163202.3 
# SBC:     163228.4 

png("plots/estim_emiBetaPrior_ZIBBregression/diagnosticPlots_unmethPop_ZIBBregression_nullModel.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_u0)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  0.006403633 
# variance   =  0.9887249 
# coef. of skewness  =  -0.1347624 
# coef. of kurtosis  =  2.593686 
# Filliben correlation coefficient  =  0.9946032 
# ******************************************************************



##################################################
##### ZIBB regression for methylated subpop ####
##################################################

show.link("ZIBB")
mod_m <- gamlss(y ~ log(med_cov), data = sub_meth, 
                sigma.formula = ~ log(med_cov), nu.formula = ~ log(med_cov), 
                family=ZIBB(nu.link = "logit", mu.link = "logit", sigma.link = "log"),
                control = gamlss.control(n.cyc = 200))
saveRDS(mod_m, file = "code/estim_emiBetaPrior_ZIBBregression/model_meth_ZIBBregression.rds")

## diagnostics
summary(mod_m)
# ******************************************************************
#   Family:  c("ZIBB", "Zero Inflated Beta Binomial") 
# 
# Call:  gamlss(formula = y ~ log(med_cov), sigma.formula = ~log(med_cov),  
#               nu.formula = ~log(med_cov), family = ZIBB(nu.link = "logit",  
#                                                         mu.link = "logit", sigma.link = "log"), data = sub_meth,  
#               control = gamlss.control(n.cyc = 200)) 
# 
# Fitting method: RS() 
# 
# ------------------------------------------------------------------
#   Mu link function:  logit
# Mu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -2.795948   0.015697 -178.12   <2e-16 ***
#   log(med_cov)  0.086802   0.004977   17.44   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Sigma link function:  log
# Sigma Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -2.28414    0.03465 -65.917  < 2e-16 ***
#   log(med_cov) -0.02996    0.01046  -2.864  0.00418 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Nu link function:  logit 
# Nu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)    -38.34    1085.36  -0.035    0.972
# log(med_cov)     1.05     360.15   0.003    0.998
# 
# ------------------------------------------------------------------
#   No. of observations in the fit:  223098 
# Degrees of Freedom for the fit:  6
# Residual Deg. of Freedom:  223092 
# at cycle:  11 
# 
# Global Deviance:     727610.2 
# AIC:     727622.2 
# SBC:     727684.1 
# ******************************************************************

png("plots/estim_emiBetaPrior_ZIBBregression/diagnosticPlots_methPop_ZIBBregression.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_m)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  0.0007623022 
# variance   =  0.9858578 
# coef. of skewness  =  0.1091545 
# coef. of kurtosis  =  3.167087 
# Filliben correlation coefficient  =  0.9982827 
# ******************************************************************

ggplot() + geom_violin(aes(as.factor(sub_meth$med_cov), mod_m$residuals))

sub_meth_smr <- sub_meth %>% 
  mutate(nu_fit = fitted(mod_m, "nu"),
         mu_fit = fitted(mod_m, "mu"),
         sigma_fit = fitted(mod_m, "sigma")) %>%
  group_by(DataSource, CellClass, SubType, med_cov) %>% 
  summarise(weight = sum(cell_MF==1) / length(cell_MF), 
            mean = mean(1-cell_MF[cell_MF<1]), 
            var = var(1-cell_MF[cell_MF<1]),
            nu_fit = max(nu_fit),
            mu_fit = max(mu_fit),
            sigma_fit = max(sigma_fit))

## MoM estimates (plug in mean and var) vs. regression estimates
sub_meth_smr %>% ggplot(aes(weight, nu_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0, 0.71) + ylim(0, 0.71) + xlab("MoM estimate of nu from IB") + ylab("Estimate of nu from ZIBB regression")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_methPop_ZIBBreg_vs_IBMoM_paramNu.png", width = 5, height = 5)
sub_meth_smr %>% ggplot(aes(1-mean, 1-mu_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.8, 1) + ylim(0.8, 1) + xlab("MoM estimate of mu from IB") + ylab("1-Estimate of mu from ZIBB regression")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_methPop_ZIBBreg_vs_IBMoM_paramMu.png", width = 5, height = 5)
sub_meth_smr %>% ggplot(aes(sqrt(var/mean/(1-mean)), sigma_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.089, 0.34) + ylim(0.089, 0.34) + xlab("MoM estimate of sigma from IB") + ylab("Estimate of sigma from ZIBB regression")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_methPop_ZIBBreg_vs_IBMoM_paramSigma.png", width = 5, height = 5)

## regression estimates vs. median coverage
sub_meth_smr %>% ggplot() + 
  geom_point(aes(log(med_cov), weight), color = "darkgrey") + 
  geom_point(aes(log(med_cov), nu_fit), color = "blue") + 
  xlab("Log median across-cell coevrage") + ylab("Ests of nu from ZIBB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_methPop_paramNu_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)
sub_meth_smr %>% ggplot() + 
  geom_point(aes(log(med_cov), 1-mean), color = "darkgrey") + 
  geom_point(aes(log(med_cov), 1-mu_fit), color = "blue") + 
  xlab("Log median across-cell coevrage") + ylab("Ests of mu from ZIBB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_methPop_paramMu_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)
sub_meth_smr %>% ggplot() + 
  geom_point(aes(log(med_cov), sqrt(var/mean/(1-mean))), color = "darkgrey") + 
  geom_point(aes(log(med_cov), sigma_fit), color = "blue") + 
  xlab("Log median across-cell coevrage") + ylab("Estimates of sigma from ZIBB (blue) and MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZIBBregression/point_methPop_paramSigma_vs_medCov_ZIBBreg&IBMoM.png", width = 5, height = 5)


## Inspectiing null model
mod_m0 <- gamlss(y ~ 1, data = sub_meth, 
                 sigma.formula = ~ 1, nu.formula = ~ 1, 
                 family=ZIBB,
                 control = gamlss.control(n.cyc = 200))
summary(mod_m0)
# Global Deviance:     727983.3 
# AIC:     727989.3 
# SBC:     728020.2 

png("plots/estim_emiBetaPrior_ZIBBregression/diagnosticPlots_methPop_ZIBBregression_nullModel.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_m0)
dev.off()



## !! Since the nu parameter in ZIBB for meth subpop is negligible, we just run BB regression

##############################################
##### BB regression for methylated subpop ####
##############################################

show.link("BB")
mod_m2 <- gamlss(y ~ log(med_cov), data = sub_meth, 
                sigma.formula = ~ log(med_cov),
                family=BB(mu.link = "logit", sigma.link = "log"),
                control = gamlss.control(n.cyc = 200))
summary(mod_m2)
# ******************************************************************
#   Family:  c("BB", "Beta Binomial") 
# 
# Call:  gamlss(formula = y ~ log(med_cov), sigma.formula = ~log(med_cov),  
#               family = BB(mu.link = "logit", sigma.link = "log"),      data = sub_meth, control = gamlss.control(n.cyc = 200)) 
# 
# 
# Fitting method: RS() 
# 
# ------------------------------------------------------------------
#   Mu link function:  logit
# Mu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -2.795939   0.016452 -169.95   <2e-16 ***
#   log(med_cov)  0.086804   0.005213   16.65   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Sigma link function:  log
# Sigma Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -2.28416    0.03475 -65.737  < 2e-16 ***
#   log(med_cov) -0.02995    0.01062  -2.821  0.00479 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   No. of observations in the fit:  223098 
# Degrees of Freedom for the fit:  4
# Residual Deg. of Freedom:  223094 
# at cycle:  2 
# 
# Global Deviance:     727610.2 
# AIC:     727618.2 
# SBC:     727659.5 
# ******************************************************************

## diagnostics
png("plots/estim_emiBetaPrior_ZIBBregression/diagnosticPlots_methPop_BBregression.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_m2)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  0.0003952794 
# variance   =  0.9866161 
# coef. of skewness  =  0.1168826 
# coef. of kurtosis  =  3.149445 
# Filliben correlation coefficient  =  0.9982746 
# ******************************************************************

sub_meth_smr <- sub_meth %>% 
  mutate(mu_fit = fitted(mod_m2, "mu"),
         sigma_fit = fitted(mod_m2, "sigma")) %>%
  group_by(DataSource, CellClass, SubType, med_cov) %>% 
  summarise(mean = mean(cell_MF), 
            var = var(cell_MF),
            mu_fit = max(mu_fit),
            sigma_fit = max(sigma_fit))

sub_meth_smr %>% ggplot(aes(mean, 1-mu_fit, color = med_cov)) + geom_point() + 
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
  geom_point(aes(log(med_cov), 1-mu_fit), color = "blue") + 
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
# Global Deviance:     727983.3 
# AIC:     727987.3 
# SBC:     728007.9

png("plots/estim_emiBetaPrior_ZABBregression/diagnosticPlots_methPop_BBregression_nullModel.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_m2_0)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  -0.001359746 
# variance   =  0.9828191 
# coef. of skewness  =  0.11851 
# coef. of kurtosis  =  3.201957 
# Filliben correlation coefficient  =  0.9981271 
# ******************************************************************