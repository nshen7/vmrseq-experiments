suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gamlss))

setwd(here::here())

## number of sites sampled from each subtype in each cluster (meth/unmeth)
sample_n <- 10000
set.seed(2022)

################################################################
# ==== Subsample subtypes in Liu2021 mice data ====
################################################################
wrapper1 <- function(info_row, cluster, cores = 16){
  dir <- "data/processed/summarized_liu2021/"
  file_list <- list.files(dir)
  file_name <- grep(pattern = paste0(sub(" ", "_", info_row$SubType), 
                                     "_", info_row$N, "cells"),
                    x = file_list, value = T)
  print(file_name)
  stopifnot("File does NOT exist." = file.exists(paste0(dir, file_name)))
  st <- fread(paste0(dir, file_name)) %>% 
    dplyr::mutate(cell_MF = cell_meth / cell_cov) %>%
    dplyr::filter(cell_cov >= 5) ## QC
  
  if (cluster == "u") {
    index <- sample(which(st$cell_MF < 0.5), sample_n)
    sub_dat <- data.frame(
      DataSource = "Liu2021",
      CellClass = info_row$CellClass,
      SubType = info_row$SubType,
      N_cell = info_row$N,
      med_cov = median(st$cell_cov),
      cell_cov = st$cell_cov[index],
      cell_meth = st$cell_meth[index],
      cell_MF = st$cell_MF[index]
    ) 
    
    ## fit ZIBB distribution on individual subtypes
    sub_dat_glm <- sub_dat %>%
      mutate(y = as.matrix(data.frame(cell_meth, cell_cov-cell_meth)))
    mod <- gamlss(y ~ 1, data = sub_dat_glm, family = ZIBB, trace = F, n.cyc = 100)
    
    sub_dat$mle_nu <- fitted(mod, "nu")[1]
    sub_dat$mle_mu <- fitted(mod, "mu")[1]
    sub_dat$mle_sigma <- unname(fitted(mod, "sigma")[1])
    
  } else if (cluster == "m"){
    index <- sample(which(st$cell_MF > 0.5), sample_n)
    sub_dat <- data.frame(
      DataSource = "Liu2021",
      CellClass = info_row$CellClass,
      SubType = info_row$SubType,
      N_cell = info_row$N,
      med_cov = median(st$cell_cov),
      cell_cov = st$cell_cov[index],
      cell_meth = st$cell_meth[index],
      cell_MF = st$cell_MF[index]
    ) 
    
    ## fit ZIBB distribution on individual subtypes
    sub_dat_glm <- sub_dat %>%
      mutate(y = as.matrix(data.frame(cell_cov-cell_meth, cell_meth)))
    mod <- gamlss(y ~ 1, data = sub_dat_glm, family = ZIBB, trace = F, n.cyc = 100)
    sub_dat$mle_nu <- fitted(mod, "nu")[1]
    sub_dat$mle_mu <- fitted(mod, "mu")[1]
    sub_dat$mle_sigma <- unname(fitted(mod, "sigma")[1])
    
    ## fit BB distribution on individual subtypes
    sub_dat_glm2 <- sub_dat %>%
      mutate(y = as.matrix(data.frame(cell_meth, cell_cov-cell_meth)))
    mod2 <- gamlss(y ~ 1, data = sub_dat_glm2, family = BB, trace = F, n.cyc = 100)
    sub_dat$bbmle_mu <- fitted(mod2, "mu")[1]
    sub_dat$bbmle_sigma <- unname(fitted(mod2, "sigma")[1])
    
  } else {
    stop("Cluster should be either 'u' or 'm'.")
  }
  
  return(sub_dat)
}


metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv") %>%
  filter(!is.na(GEO_accession) & !is.na(SubType) & !is.na(FilePath))
subtype_smr <- metadata[, .(.N), by = .(CellClass, SubType)] %>% 
  arrange(desc(N)) %>% 
  filter(!grepl("Outlier", SubType)) %>%
  mutate(Bin = cut(N, breaks = seq(0, ceiling(max(N)/100)*100, 100), labels = F)) %>%
  group_by(Bin) %>% mutate(nSubTypeInBin = n())

## Select (at most) 3 subtypes from each bin of 100 in N cell
set.seed(2022)
subtype_smr_sel <- rbind(subtype_smr %>% 
                           filter(nSubTypeInBin <= 3),
                         subtype_smr %>% 
                           filter(nSubTypeInBin > 3) %>%
                           group_by(Bin) %>% 
                           sample_n(3) %>% 
                           mutate(nSubTypeInBin = n())
) %>% 
  arrange(desc(N))

# info <- subtype_smr_sel
info <- subtype_smr_sel %>% filter(N >= 100)

## compute density and fit parameters
sub_liu2021_u <- do.call(rbind, map(1:nrow(info), ~ wrapper1(info[.x,], 'u')))
fwrite(sub_liu2021_u, "data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_liu2021_clusterUnmeth.csv")
sub_liu2021_m <- do.call(rbind, map(1:nrow(info), ~ wrapper1(info[.x,], 'm')))
fwrite(sub_liu2021_m, "data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_liu2021_clusterMeth.csv")


################################################################
# ==== Subsample subtypes in Luo2017 mice data ====
################################################################

wrapper2 <- function(info_row, cluster, cores = 16){
  dir <- "data/processed/processed_luo2017_mice/"
  file_list <- list.files(dir)
  file_name <- grep(pattern = paste0("subtype_", sub("/", "", info_row$SubType), 
                                     ".*_qced$"),
                    x = file_list, value = T)
  print(file_name)
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))
  
  if (cluster == "u") {
    index <- sample(which(values(cells.se)$cell_MF < 0.5), sample_n)
    sub_dat <- data.frame(
      DataSource = "Luo2017m",
      CellClass = info_row$CellClass,
      SubType = info_row$SubType,
      N_cell = info_row$N,
      med_cov = median(values(cells.se)$cell_cov),
      cell_cov = values(cells.se)$cell_cov[index],
      cell_meth = values(cells.se)$cell_meth[index],
      cell_MF = values(cells.se)$cell_MF[index]
    ) 
    
    ## fit ZIBB distribution on individual subtypes
    sub_dat_glm <- sub_dat %>%
      mutate(y = as.matrix(data.frame(cell_meth, cell_cov-cell_meth)))
    mod <- gamlss(y ~ 1, data = sub_dat_glm, family = ZIBB, trace = F, n.cyc = 100)
    
    sub_dat$mle_nu <- fitted(mod, "nu")[1]
    sub_dat$mle_mu <- fitted(mod, "mu")[1]
    sub_dat$mle_sigma <- unname(fitted(mod, "sigma")[1])
    
  } else if (cluster == "m"){
    index <- sample(which(values(cells.se)$cell_MF > 0.5), sample_n)
    sub_dat <- data.frame(
      DataSource = "Luo2017m",
      CellClass = info_row$CellClass,
      SubType = info_row$SubType,
      N_cell = info_row$N,
      med_cov = median(values(cells.se)$cell_cov),
      cell_cov = values(cells.se)$cell_cov[index],
      cell_meth = values(cells.se)$cell_meth[index],
      cell_MF = values(cells.se)$cell_MF[index]
    ) 
    
    ## fit ZIBB distribution on individual subtypes
    sub_dat_glm <- sub_dat %>%
      mutate(y = as.matrix(data.frame(cell_cov-cell_meth, cell_meth)))
    mod <- gamlss(y ~ 1, data = sub_dat_glm, family = ZIBB, trace = F, n.cyc = 100)
    sub_dat$mle_nu <- fitted(mod, "nu")[1]
    sub_dat$mle_mu <- fitted(mod, "mu")[1]
    sub_dat$mle_sigma <- unname(fitted(mod, "sigma")[1])
    
    ## fit BB distribution on individual subtypes
    sub_dat_glm2 <- sub_dat %>%
      mutate(y = as.matrix(data.frame(cell_meth, cell_cov-cell_meth)))
    mod2 <- gamlss(y ~ 1, data = sub_dat_glm2, family = BB, trace = F, n.cyc = 100)
    sub_dat$bbmle_mu <- fitted(mod2, "mu")[1]
    sub_dat$bbmle_sigma <- unname(fitted(mod2, "sigma")[1])
    
  } else {
    stop("Cluster should be either 'u' or 'm'.")
  }
  
  return(sub_dat)
}

## information of all subtypes used for training
metadata <- fread("../../DXM_extend_chr1/data/metadata/sample_info_processed.csv")
subtype_smr <- metadata[, .(.N), by = .(Neuron_type1, Neuron_type3)] %>%
  arrange(desc(N)) %>%
  filter(N >= 100) %>%
  dplyr::rename(CellClass = Neuron_type1, SubType = Neuron_type3) %>%
  mutate(CellClass = recode(CellClass, "Excitatory" = "Exc", "Inhibitory" = "Inh")) %>%
  mutate(Bin = cut(N, breaks = seq(0, ceiling(max(N)/100)*100, 100), labels = F)) %>%
  group_by(Bin) %>% mutate(nSubTypeInBin = n())

## Select (at most) 3 subtypes from each bin of 100 in N cell
set.seed(2023)
subtype_smr_sel <- rbind(subtype_smr %>% 
                           filter(nSubTypeInBin <= 3),
                         subtype_smr %>% 
                           filter(nSubTypeInBin > 3) %>%
                           group_by(Bin) %>% 
                           sample_n(3) %>% 
                           mutate(nSubTypeInBin = n())
) %>% 
  arrange(desc(N))

info <- subtype_smr_sel

## compute density and fit parameters
sub_luo2017m_u <- do.call(rbind, map(1:nrow(info), ~wrapper2(info[.x,], 'u')))
fwrite(sub_luo2017m_u, "data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017mice_clusterUnmeth.csv")
sub_luo2017m_m <- do.call(rbind, map(1:nrow(info), ~wrapper2(info[.x,], 'm')))
fwrite(sub_luo2017m_m, "data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017mice_clusterMeth.csv")


################################################################
# ==== Subsample subtypes in Luo2017 human data ====
################################################################

wrapper3 <- function(info_row, cluster, cores = 16){
  dir <- "data/processed/processed_luo2017_human/"
  file_list <- list.files(dir)
  file_name <- grep(pattern = paste0("subtype_", sub("/", "", info_row$SubType), 
                                     ".*_qced$"),
                    x = file_list, value = T)
  print(file_name)
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))

  if (cluster == "u") {
    index <- sample(which(values(cells.se)$cell_MF < 0.5), sample_n)
    sub_dat <- data.frame(
      DataSource = "Luo2017h",
      CellClass = NA,
      SubType = info_row$SubType,
      N_cell = info_row$N,
      med_cov = median(values(cells.se)$cell_cov),
      cell_cov = values(cells.se)$cell_cov[index],
      cell_meth = values(cells.se)$cell_meth[index],
      cell_MF = values(cells.se)$cell_MF[index]
    ) 
    
    ## fit ZIBB distribution on individual subtypes
    sub_dat_glm <- sub_dat %>%
      dplyr::select(-CellClass) %>%
      mutate(y = as.matrix(data.frame(cell_meth, cell_cov-cell_meth)))
    mod <- gamlss(y ~ 1, data = sub_dat_glm, family = ZIBB, trace = F, n.cyc = 100)
    
    sub_dat$mle_nu <- fitted(mod, "nu")[1]
    sub_dat$mle_mu <- fitted(mod, "mu")[1]
    sub_dat$mle_sigma <- unname(fitted(mod, "sigma")[1])
    
  } else if (cluster == "m"){
    index <- sample(which(values(cells.se)$cell_MF > 0.5), sample_n)
    sub_dat <- data.frame(
      DataSource = "Luo2017h",
      CellClass = NA,
      SubType = info_row$SubType,
      N_cell = info_row$N,
      med_cov = median(values(cells.se)$cell_cov),
      cell_cov = values(cells.se)$cell_cov[index],
      cell_meth = values(cells.se)$cell_meth[index],
      cell_MF = values(cells.se)$cell_MF[index]
    ) 
    
    ## fit ZIBB distribution on individual subtypes
    sub_dat_glm <- sub_dat %>%
      dplyr::select(-CellClass) %>%
      mutate(y = as.matrix(data.frame(cell_cov-cell_meth, cell_meth)))
    mod <- gamlss(y ~ 1, data = sub_dat_glm, family = ZIBB, trace = T, n.cyc = 200)
    sub_dat$mle_nu <- fitted(mod, "nu")[1]
    sub_dat$mle_mu <- fitted(mod, "mu")[1]
    sub_dat$mle_sigma <- unname(fitted(mod, "sigma")[1])
    
    ## fit BB distribution on individual subtypes
    sub_dat_glm2 <- sub_dat %>%
      dplyr::select(-CellClass) %>%
      mutate(y = as.matrix(data.frame(cell_meth, cell_cov-cell_meth)))
    mod2 <- gamlss(y ~ 1, data = sub_dat_glm2, family = BB, trace = T, n.cyc = 200)
    sub_dat$bbmle_mu <- fitted(mod2, "mu")[1]
    sub_dat$bbmle_sigma <- unname(fitted(mod2, "sigma")[1])
   
  } else {
    stop("Cluster should be either 'u' or 'm'.")
  }
  
  return(sub_dat)
}

## information of all subtypes used for training
metadata <- fread("data/metadata/metadata_luo2017/NIHMS893063-supplement-Table_S2_csv.csv", skip = 1)[, .(Sample, `Neuron type`)]
subtype_smr <- metadata[, .(.N), by = .(`Neuron type`)] %>%
  arrange(desc(N)) %>%
  filter(N >= 100) %>%
  dplyr::rename(SubType = `Neuron type`) %>%
  mutate(Bin = cut(N, breaks = seq(0, ceiling(max(N)/100)*100, 100), labels = F)) %>%
  group_by(Bin) %>% mutate(nSubTypeInBin = n())

## Select (at most) 3 subtypes from each bin of 100 in N cell
set.seed(2024)
subtype_smr_sel <- rbind(subtype_smr %>% 
                           filter(nSubTypeInBin <= 3),
                         subtype_smr %>% 
                           filter(nSubTypeInBin > 3) %>%
                           group_by(Bin) %>% 
                           sample_n(3) %>% 
                           mutate(nSubTypeInBin = n())
) %>% 
  arrange(desc(N))

info <- subtype_smr_sel

## compute density and fit parameters
sub_luo2017h_u <- do.call(rbind, map(1:nrow(info), ~wrapper3(info[.x,], 'u')))
fwrite(sub_luo2017h_u, "data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017human_clusterUnmeth.csv")
sub_luo2017h_m <- do.call(rbind, map(1:nrow(info), ~wrapper3(info[.x,], 'm')))
fwrite(sub_luo2017h_m, "data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017human_clusterMeth.csv")



##############
##### EDA ####
##############

## Unmethylated cluster
sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_liu2021_clusterUnmeth.csv")
sub_luo2017m <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017mice_clusterUnmeth.csv")
sub_luo2017h <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017human_clusterUnmeth.csv")
sub_unmeth_smr <- rbind(sub_liu2021, sub_luo2017m, sub_luo2017h) %>% 
  group_by(DataSource, CellClass, SubType, med_cov) %>%
  summarise(nu_mle = max(mle_nu),
            mu_mle = max(mle_mu),
            sigma_mle = max(mle_sigma))
sub_unmeth_smr$CellClass[is.na(sub_unmeth_smr$CellClass)] <- "Unknown"

sub_unmeth_smr %>%
  ggplot(aes(log(med_cov), nu_mle)) +
  geom_smooth(method = "loess", method.args = list(degree = 1)) +
  geom_point(aes(color = CellClass, shape = DataSource)) +
  scale_shape_manual(values=c(16, 2, 3))
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_MLE_vs_medCov_unmeth_paramNu.png", width = 7, height = 5)
sub_unmeth_smr %>%
  ggplot(aes(log(med_cov), mu_mle)) +
  geom_smooth(method = "loess", method.args = list(degree = 1)) +
  geom_point(aes(color = CellClass, shape = DataSource)) +
  scale_shape_manual(values=c(16, 2, 3))
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_MLE_vs_medCov_unmeth_paramMu.png", width = 7, height = 5)
sub_unmeth_smr %>%
  ggplot(aes(log(med_cov), sigma_mle)) +
  geom_smooth(method = "loess", method.args = list(degree = 1)) +
  geom_point(aes(color = CellClass, shape = DataSource)) +
  scale_shape_manual(values=c(16, 2, 3))
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_MLE_vs_medCov_unmeth_paramSigma.png", width = 7, height = 5)


## Methylated cluster
sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_liu2021_clusterMeth.csv")
sub_luo2017m <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017mice_clusterMeth.csv")
sub_luo2017h <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017human_clusterMeth.csv")
sub_meth_smr <- rbind(sub_liu2021, sub_luo2017m, sub_luo2017h) %>% 
  group_by(DataSource, CellClass, SubType, med_cov) %>%
  summarise(nu_mle = max(mle_nu),
            mu_mle = max(mle_mu),
            sigma_mle = max(mle_sigma),
            mu_bbmle = max(bbmle_mu),
            sigma_bbmle = max(bbmle_sigma)
            )
sub_meth_smr$CellClass[is.na(sub_meth_smr$CellClass)] <- "Unknown"

sub_meth_smr %>%
  ggplot(aes(log(med_cov), nu_mle)) +
  geom_smooth(method = "loess", method.args = list(degree = 1)) +
  geom_point(aes(color = CellClass, shape = DataSource)) +
  scale_shape_manual(values=c(16, 2, 3))
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_MLE_vs_medCov_meth_paramNu.png", width = 7, height = 5)
sub_meth_smr %>%
  ggplot(aes(log(med_cov), mu_mle)) +
  geom_smooth(method = "loess", method.args = list(degree = 1)) +
  geom_point(aes(color = CellClass, shape = DataSource)) +
  scale_shape_manual(values=c(16, 2, 3))
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_MLE_vs_medCov_meth_paramMu.png", width = 7, height = 5)
sub_meth_smr %>%
  ggplot(aes(log(med_cov), sigma_mle)) +
  geom_smooth(method = "loess", method.args = list(degree = 1)) +
  geom_point(aes(color = CellClass, shape = DataSource)) +
  scale_shape_manual(values=c(16, 2, 3))
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_MLE_vs_medCov_meth_paramSigma.png", width = 7, height = 5)


sub_meth_smr %>%
  ggplot(aes(log(med_cov), mu_bbmle)) +
  geom_smooth(method = "loess", method.args = list(degree = 1)) +
  geom_point(aes(color = CellClass, shape = DataSource)) +
  scale_shape_manual(values=c(16, 2, 3))
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_BBMLE_vs_medCov_meth_paramMu.png", width = 7, height = 5)
sub_meth_smr %>%
  ggplot(aes(log(med_cov), sigma_bbmle)) +
  geom_smooth(method = "loess", method.args = list(degree = 1)) +
  geom_point(aes(color = CellClass, shape = DataSource)) +
  scale_shape_manual(values=c(16, 2, 3))
ggsave("plots/estim_emiBetaPrior_ZIBBregression2/point_BBMLE_vs_medCov_meth_paramSigma.png", width = 7, height = 5)
