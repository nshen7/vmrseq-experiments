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
sub_meth$y <- sub_meth$y[,2:1]

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
# (Intercept)  -2.805258   0.007357 -381.33   <2e-16 ***
#   log(med_cov)  0.088394   0.002331   37.92   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Sigma link function:  log
# Sigma Coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)  -2.303546   0.015620 -147.470  < 2e-16 ***
#   log(med_cov) -0.026199   0.004774   -5.488 4.07e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# ------------------------------------------------------------------
#   Nu link function:  logit 
# Nu Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)    -22.23     131.07  -0.170    0.865
# log(med_cov)    -3.66     131.07  -0.028    0.978
# 
# ------------------------------------------------------------------
#   No. of observations in the fit:  1116502 
# Degrees of Freedom for the fit:  6
# Residual Deg. of Freedom:  1116496 
# at cycle:  12 
# 
# Global Deviance:     3634069 
# AIC:     3634081 
# SBC:     3634152 
# ******************************************************************
#   
png("plots/estim_emiBetaPrior_ZIBBregression/diagnosticPlots_methPop_ZIBBregression.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_m)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  0.0006295749 
# variance   =  0.9838288 
# coef. of skewness  =  0.1152611 
# coef. of kurtosis  =  3.150013 
# Filliben correlation coefficient  =  0.9983521 
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
# Global Deviance:     3635965 
# AIC:     3635971 
# SBC:     3636007 

png("plots/estim_emiBetaPrior_ZIBBregression/diagnosticPlots_methPop_ZIBBregression_nullModel.png", width = 1000, height = 700)
par(mar=c(2,1.5,1.5,1.5)); plot(mod_m0)
dev.off()
# ******************************************************************
#   Summary of the Randomised Quantile Residuals
# mean   =  -0.00286046 
# variance   =  0.9836951 
# coef. of skewness  =  0.1197339 
# coef. of kurtosis  =  3.188275 
# Filliben correlation coefficient  =  0.9982283 
# ******************************************************************


