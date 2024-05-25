source("code/SETPATHS.R")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BiocParallel))
devtools::load_all("../vmrseq-package/vmrseq/")
set.seed(2022)

###################################################
######## Train emission beta-binomial model #######
###################################################
dir <- "data/interim/estim_emiBetaPrior_ZIBBregression2/"
sub_liu2021 <- fread(paste0(dir, "emiBetaPrior_subtype_subsample_liu2021_clusterMeth.csv"))
# sub_luo2017m <- fread(paste0(dir, "emiBetaPrior_subtype_subsample_luo2017mice_clusterMeth.csv"))
# sub_luo2017h <- fread(paste0(dir, "emiBetaPrior_subtype_subsample_luo2017human_clusterMeth.csv"))

sub_meth <- sub_liu2021 %>%
  as_tibble() %>%
  mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth))) %>%
  dplyr::select(-CellClass) # CellClass contains NA from liu2017 human data so removed

mod_m <- gamlss(y ~ 1, data = sub_meth,
                sigma.formula = ~ 1,
                family = BB(mu.link = "logit", sigma.link = "log"),
                n.cyc = 50)
params_m <- data.frame(
  mu = fitted(mod_m, "mu")[1],
  sigma = fitted(mod_m, "sigma")[1]
)

sub_meth_smr <- sub_meth %>%
  mutate(mu_fit = fitted(mod_m, "mu"),
         sigma_fit = fitted(mod_m, "sigma")) %>%
  group_by(DataSource, SubType, med_cov) %>%
  summarise(mu_mle = max(bbmle_mu),
            sigma_mle = max(bbmle_sigma),
            mu_fit = max(mu_fit),
            sigma_fit = max(sigma_fit))
fwrite(sub_meth_smr, "code/params_sens_analysis/fitted_params_beta_prior_m.csv")

## regression estimates vs. median coverage
sub_meth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), mu_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), mu_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Ests of mu from BB (blue) & MoM (grey)") +
  theme_classic()
ggsave("code/params_sens_analysis/point_methPop_paramMu_vs_medCov_BBreg&mle.png", width = 5, height = 5)
sub_meth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), sigma_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), sigma_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Estimates of sigma from BB (blue) and MoM (grey)") +
  theme_classic()
ggsave("code/params_sens_analysis/point_methPop_paramSigma_vs_medCov_BBreg&mle.png", width = 5, height = 5)



sub_liu2021 <- fread(paste0(dir, "emiBetaPrior_subtype_subsample_liu2021_clusterUnmeth.csv"))
# sub_luo2017m <- fread(paste0(dir, "emiBetaPrior_subtype_subsample_luo2017mice_clusterUnmeth.csv"))
# sub_luo2017h <- fread(paste0(dir, "emiBetaPrior_subtype_subsample_luo2017human_clusterUnmeth.csv"))
sub_unmeth <- sub_liu2021 %>%
  as_tibble() %>%
  mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth))) %>%
  dplyr::select(-CellClass) # CellClass contains NA from liu2017 human data so removed

## Fit ZIBB regression on nu, mu, sigma against log(med_cov)
show.link("ZIBB")
mod_u <- gamlss(y ~ cs(log(med_cov)), data = sub_unmeth,
                sigma.formula = ~ cs(log(med_cov)),
                nu.formula = ~ cs(log(med_cov)),
                family = ZIBB(nu.link = "logit", mu.link = "logit", sigma.link = "log"),
                n.cyc = 100)
med_cov <- 1:400
params_u <- data.frame(
  med_cov = med_cov,
  nu = predict(mod_u, what = "nu", type = "response", newdata = data.frame(med_cov = med_cov)),
  mu = predict(mod_u, what = "mu", type = "response", newdata = data.frame(med_cov = med_cov)),
  sigma = predict(mod_u, what = "sigma", type = "response", newdata = data.frame(med_cov = med_cov))
)

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
fwrite(sub_unmeth_smr, "code/params_sens_analysis/fitted_params_beta_prior_u.csv")

## regression estimates vs. median coverage
sub_unmeth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), nu_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), nu_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Ests of nu from ZIBB (blue) & MoM (grey)") +
  theme_classic()
ggsave("code/params_sens_analysis/point_unmeth_paramNu_vs_medCov_ZIBBreg&mle.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), mu_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), mu_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Ests of mu from ZIBB (blue) & MoM (grey)") +
  theme_classic()
ggsave("code/params_sens_analysis/point_unmeth_paramMu_vs_medCov_ZIBBreg&mle.png", width = 5, height = 5)
sub_unmeth_smr %>% ggplot() +
  geom_point(aes(log(med_cov), sigma_mle), color = "darkgrey") +
  geom_point(aes(log(med_cov), sigma_fit), color = "blue") +
  xlab("Log median across-cell coevrage") + ylab("Estimates of sigma from ZIBB (blue) and MoM (grey)") +
  theme_classic()
ggsave("code/params_sens_analysis/point_unmeth_paramSigma_vs_medCov_ZIBBreg&mle.png", width = 5, height = 5)


########################################################################
######## Train transition probablities used as default in vmrseq #######
########################################################################

# Using default parameters for estimating transition probs:
max_dist_bp = 2000; buffer_bp = 3000
lags = 1:10
degree = 2; span = 0.02

##### Summarize info in Liu2021 subtypes and form into list of data.frames

wrapper1 <- function(subtype_dir){

  st <- fread(subtype_dir) %>%
    filter(cell_cov >= 5) %>% # QC: across-cell coverage >= 5
    mutate(cell_MF = cell_meth / cell_cov) %>%
    mutate(state = round(cell_MF)) %>%
    group_by(chr) %>%
    mutate(lag_state = lag(state, 1))

  # decide the states for sites with MF = 0.5 based on its previous site
  st$state[which(st$cell_MF==0.5)] <- st$lag_state[which(st$cell_MF==0.5)]

  smr <- vmrseq:::.computeProb1Unit(
    df = st %>% select(chr, pos, state),
    max_dist_bp = max_dist_bp,
    buffer_bp = buffer_bp,
    lags = lags
  )

  return(smr)
}

# Information of all subtypes used for training
dir <- ("data/processed/summarized_liu2021/")
subtype_dirs <- paste0(dir, list.files(dir))

# Compute empirical probabilities from each subtype
smr_liu <- do.call(
  rbind,
  bplapply(subtype_dirs, function(.x) wrapper1(.x))
)
print("Finished processing Liu2021 data.")

smr <- smr_liu

tp0 <- vmrseq:::.estimTransitProbsFromSummary(
  smr_units = smr,
  max_dist_bp = max_dist_bp,
  buffer_bp = buffer_bp,
  degree = degree,
  span = span
)
saveRDS(tp0, "code/params_sens_analysis/tp0.rds")

# plot loess-fitted transition probs
tp.plot(tp0)
ggsave("code/params_sens_analysis/default_transit_probs.png", width = 8, height = 6)

