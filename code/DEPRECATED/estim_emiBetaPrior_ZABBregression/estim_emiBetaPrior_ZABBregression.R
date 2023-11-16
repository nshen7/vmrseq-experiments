library(tidyverse)
library(gamlss)
library(gamlss.tr)
setwd(here::here())

sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_liu2021.csv")
sub_luo2017 <- fread("data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_luo2017.csv")
sub_unmeth <- rbind(sub_luo2017, sub_liu2021) %>% 
  dplyr::filter(cell_MF < 0.5) %>%
  mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth)))
sub_meth <- rbind(sub_luo2017, sub_liu2021) %>% 
  dplyr::filter(cell_MF > 0.5) %>%
  mutate(y = as.matrix(data.frame(cell_cov - cell_meth, cell_meth)))

########################################
##### regression for unmethylated subpop

show.link("ZABB")
mod_u <- gamlss(y ~ log(med_cov), data = sub_unmeth, 
                sigma.formula = ~ log(med_cov), nu.formula = ~ log(med_cov), 
                family=ZABB(nu.link = "logit", mu.link = "logit", sigma.link = "log"),
                control = gamlss.control(n.cyc = 200))
## diagnostics
summary(mod_u)
# ******************************************************************
#   Family:  c("ZABB", "Zero Adjusted Beta Binomial") 
# 
# Call:  gamlss(formula = y ~ log(med_cov), sigma.formula = ~log(med_cov),  
#               nu.formula = ~log(med_cov), family = ZABB(nu.link = "logit",  
#                                                         mu.link = "logit", sigma.link = "log"), data = sub_unmeth,  
#               control = gamlss.control(n.cyc = 200, mu.step = 1,  
#                                        sigma.step = 1)) 
# 
# Fitting method: RS() 
# 
# ------------------------------------------------------------------
#   Mu link function:  logit
# Mu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.51657    0.06379  -8.097 5.75e-16 ***
#   log(med_cov) -0.43136    0.02220 -19.432  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Sigma link function:  log
# Sigma Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -4.01762    0.11381  -35.30   <2e-16 ***
#   log(med_cov)  0.76847    0.03419   22.47   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Nu link function:  logit 
# Nu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.40507    0.04797   29.29   <2e-16 ***
#   log(med_cov) -0.44704    0.01577  -28.35   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   No. of observations in the fit:  45370 
# Degrees of Freedom for the fit:  6
# Residual Deg. of Freedom:  45364 
# at cycle:  7 
# 
# Global Deviance:     161746 
# AIC:     161758 
# SBC:     161810.3 
# ******************************************************************
  
png("plots/estim_emiBetaPrior_ZABBregression/diagnosticPlots_unmethPop_ZABBregression.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_u)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  0.0006083936 
# variance   =  0.9911001 
# coef. of skewness  =  -0.1171253 
# coef. of kurtosis  =  2.59584 
# Filliben correlation coefficient  =  0.9945059 
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
  xlim(0.3, 0.71) + ylim(0.3, 0.71) + xlab("MoM estimate of nu from IB") + ylab("Estimate of nu from ZABB regression")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_unmethPop_ZABBreg_vs_IBMoM_paramNu.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot(aes(mean, mu_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.08, 0.26) + ylim(0.08, 0.26) + xlab("MoM estimate of mu from IB") + ylab("Estimate of mu from ZABB regression")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_unmethPop_ZABBreg_vs_IBMoM_paramMu.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot(aes(sqrt(var/mean/(1-mean)), sigma_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.05, 0.5) + ylim(0.05, 0.5) + xlab("MoM estimate of sigma from IB") + ylab("Estimate of sigma from ZABB regression")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_unmethPop_ZABBreg_vs_IBMoM_paramSigma.png", width = 5, height = 5)

## regression estimates vs. median coverage
sub_unmeth_smr %>% ggplot() + 
  geom_point(aes(log(med_cov), weight), color = "darkgrey") + 
  geom_point(aes(log(med_cov), nu_fit), color = "blue") + 
  xlab("Log median across-cell coevrage") + ylab("Ests of nu from ZABB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_unmethPop_paramNu_vs_medCov_ZABBreg&IBMoM.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot() + 
  geom_point(aes(log(med_cov), mean), color = "darkgrey") + 
  geom_point(aes(log(med_cov), mu_fit), color = "blue") + 
  xlab("Log median across-cell coevrage") + ylab("Ests of mu from ZABB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_unmethPop_paramMu_vs_medCov_ZABBreg&IBMoM.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot() + 
  geom_point(aes(log(med_cov), sqrt(var/mean/(1-mean))), color = "darkgrey") + 
  geom_point(aes(log(med_cov), sigma_fit), color = "blue") + 
  xlab("Log median across-cell coevrage") + ylab("Estimates of sigma from ZABB (blue) and MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_unmethPop_paramSigma_vs_medCov_ZABBreg&IBMoM.png", width = 5, height = 5)

## Inspecting null model
mod_u0 <- gamlss(y ~ 1, data = sub_unmeth, 
                 sigma.formula = ~ 1, nu.formula = ~ 1, 
                 family=ZABB,
                 control = gamlss.control(n.cyc = 200))
summary(mod_u0)
# Global Deviance:     163196.3 
# AIC:     163202.3 
# SBC:     163228.4 

png("plots/estim_emiBetaPrior_ZABBregression/diagnosticPlots_unmethPop_ZABBregression_nullModel.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_u0)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  0.008713003 
# variance   =  0.9834611 
# coef. of skewness  =  -0.1271453 
# coef. of kurtosis  =  2.587371 
# Filliben correlation coefficient  =  0.9945435 
# ******************************************************************

########################################
##### regression for methylated subpop
show.link("ZABB")
mod_m <- gamlss(y ~ log(med_cov), data = sub_meth, 
                sigma.formula = ~ log(med_cov), nu.formula = ~ log(med_cov), 
                family=ZABB(nu.link = "logit", mu.link = "logit", sigma.link = "log"),
                control = gamlss.control(n.cyc = 200))
## diagnostics
summary(mod_m)
# ******************************************************************
#   Family:  c("ZABB", "Zero Adjusted Beta Binomial") 
# 
# Call:  gamlss(formula = y ~ log(med_cov), sigma.formula = ~log(med_cov),  
#               nu.formula = ~log(med_cov), family = ZABB(nu.link = "logit",  
#                                                         mu.link = "logit", sigma.link = "log"), data = sub_meth,  
#               control = gamlss.control(n.cyc = 200)) 
# 
# Fitting method: RS() 
# 
# ------------------------------------------------------------------
#   Mu link function:  logit
# Mu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -4.87732    0.20736 -23.521   <2e-16 ***
#   log(med_cov)  0.48896    0.05654   8.649   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Sigma link function:  log
# Sigma Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -1.76357    0.06098 -28.921   <2e-16 ***
#   log(med_cov) -0.03407    0.01792  -1.901   0.0573 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Nu link function:  logit 
# Nu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.91082    0.02375   122.6   <2e-16 ***
#   log(med_cov) -1.09390    0.00803  -136.2   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   No. of observations in the fit:  223098 
# Degrees of Freedom for the fit:  6
# Residual Deg. of Freedom:  223092 
# at cycle:  24 
# 
# Global Deviance:     729952.9 
# AIC:     729964.9 
# SBC:     730026.8 
# ******************************************************************

png("plots/estim_emiBetaPrior_ZABBregression/diagnosticPlots_methPop_ZABBregression.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_m)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  0.00144442 
# variance   =  0.9949318 
# coef. of skewness  =  -0.008355331 
# coef. of kurtosis  =  2.897266 
# Filliben correlation coefficient  =  0.9992767 
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
  xlim(0.1, 0.71) + ylim(0.1, 0.71) + xlab("MoM estimate of nu from IB") + ylab("Estimate of nu from ZABB regression")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_methPop_ZABBreg_vs_IBMoM_paramNu.png", width = 5, height = 5)
sub_meth_smr %>% ggplot(aes(1-mean, 1-mu_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.8, 1) + ylim(0.8, 1) + xlab("MoM estimate of mu from IB") + ylab("1-Estimate of mu from ZABB regression")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_methPop_ZABBreg_vs_IBMoM_paramMu.png", width = 5, height = 5)
sub_meth_smr %>% ggplot(aes(sqrt(var/mean/(1-mean)), sigma_fit, color = med_cov)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlim(0.147, 0.34) + ylim(0.147, 0.34) + xlab("MoM estimate of sigma from IB") + ylab("Estimate of sigma from ZABB regression")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_methPop_ZABBreg_vs_IBMoM_paramSigma.png", width = 5, height = 5)

## regression estimates vs. median coverage
sub_meth_smr %>% ggplot() + 
  geom_point(aes(log(med_cov), weight), color = "darkgrey") + 
  geom_point(aes(log(med_cov), nu_fit), color = "blue") + 
  xlab("Log median across-cell coevrage") + ylab("Ests of nu from ZABB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_methPop_paramNu_vs_medCov_ZABBreg&IBMoM.png", width = 5, height = 5)
sub_meth_smr %>% ggplot() + 
  geom_point(aes(log(med_cov), 1-mean), color = "darkgrey") + 
  geom_point(aes(log(med_cov), 1-mu_fit), color = "blue") + 
  xlab("Log median across-cell coevrage") + ylab("Ests of mu from ZABB (blue) & MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_methPop_paramMu_vs_medCov_ZABBreg&IBMoM.png", width = 5, height = 5)
sub_meth_smr %>% ggplot() + 
  geom_point(aes(log(med_cov), sqrt(var/mean/(1-mean))), color = "darkgrey") + 
  geom_point(aes(log(med_cov), sigma_fit), color = "blue") + 
  xlab("Log median across-cell coevrage") + ylab("Estimates of sigma from ZABB (blue) and MoM (grey)")
ggsave("plots/estim_emiBetaPrior_ZABBregression/point_methPop_paramSigma_vs_medCov_ZABBreg&IBMoM.png", width = 5, height = 5)


## Inspectiing null model
mod_m0 <- gamlss(y ~ 1, data = sub_meth, 
                 sigma.formula = ~ 1, nu.formula = ~ 1, 
                 family=ZABB,
                 control = gamlss.control(n.cyc = 200))
summary(mod_m0)
# Global Deviance:     750593.9 
# AIC:     750599.9 
# SBC:     750630.9 

png("plots/estim_emiBetaPrior_ZABBregression/diagnosticPlots_methPop_ZABBregression_nullModel.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_m0)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  -0.001771871 
# variance   =  0.9984826 
# coef. of skewness  =  -0.01592487 
# coef. of kurtosis  =  2.942641 
# Filliben correlation coefficient  =  0.9995339 
# ******************************************************************