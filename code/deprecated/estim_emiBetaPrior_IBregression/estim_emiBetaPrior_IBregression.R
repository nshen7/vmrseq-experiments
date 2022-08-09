library(tidyverse)
library(gamlss)
setwd(here::here())

## number of sites sampled from each subtype
sample_n <- 10000
set.seed(2022)

################################################################
# ==== Compute densities and fit parameters in Liu2021 data ====
################################################################
wrapper1 <- function(info_row, cores = 16){
  dir <- "data/processed/processed_liu2021/"
  file_list <- list.files(dir)
  file_name <- grep(pattern = paste0("sample", info_row$Sample,
                                     "_", info_row$GEO_accession, 
                                     "_subtype_", sub(" ", "_", info_row$SubType),
                                     ".*_qced$"),
                    x = file_list, value = T)
  print(file_name)
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))
  
  index <- sample(1:nrow(cells.se), sample_n)
  sub_dat <- data.frame(
    DataSource = "Liu2021",
    CellClass = info_row$CellClass,
    SubType = info_row$SubType, 
    N_cell = ncol(cells.se),
    med_cov = median(values(cells.se)$cell_cov),
    cell_cov = values(cells.se)$cell_cov[index],
    cell_meth = values(cells.se)$cell_meth[index],
    cell_MF = values(cells.se)$cell_MF[index]
  )
  
  return(sub_dat)
}

## information of all subtypes used for training
metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv")
cell_counts <- metadata[, .(.N), by = .(Sample, CellClass, SubType, GEO_accession)]

cell_counts_exc <- cell_counts %>% # selected subtypes in excitatory cells for training
  filter(CellClass=="Exc") %>%
  arrange(desc(N)) %>%
  filter(!is.na(SubType), N >= 100) %>%
  dplyr::slice(c(5L, which(!duplicated(SubType)&!duplicated(Sample)))) %>%
  dplyr::slice(-2L)

cell_counts_inh <- cell_counts %>%
  filter(CellClass=="Inh") %>%
  arrange(desc(N)) %>%
  filter(!is.na(SubType), N >= 100) %>%
  dplyr::slice(which(!duplicated(SubType)&!duplicated(Sample)))

info <- rbind(cell_counts_exc,cell_counts_inh)

## compute density and fit parameters
sub_liu2021 <- do.call(rbind, map(1:nrow(info), function(i) wrapper1(info[i])))
fwrite(sub_liu2021, "data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_liu2021.csv")




################################################################
# ==== Compute densities and fit parameters in Luo2017 data ====
################################################################

wrapper2 <- function(info_row, cores = 16){
  dir <- "data/processed/processed_luo2017_mice/"
  file_list <- list.files(dir)
  file_name <- grep(pattern = paste0("subtype_", sub("/", "", info_row$SubType), ".*_qced$"),
                    x = file_list, value = T)
  print(file_name)
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))
  
  index <- sample(1:nrow(cells.se), sample_n)
  sub_dat <- data.frame(
    DataSource = "Luo2017",
    CellClass = info_row$CellClass,
    SubType = info_row$SubType, 
    N_cell = ncol(cells.se),
    med_cov = median(values(cells.se)$cell_cov),
    cell_cov = values(cells.se)$cell_cov[index],
    cell_meth = values(cells.se)$cell_meth[index],
    cell_MF = values(cells.se)$cell_MF[index]
  )
  
  return(sub_dat)
}

## information of all subtypes used for training
metadata <- fread("../../DXM_extend_chr1/data/metadata/sample_info_processed.csv")
info <- metadata[, .(.N), by = .(Neuron_type1, Neuron_type3)] %>% 
  arrange(desc(N)) %>%
  filter(N >= 100) %>%
  dplyr::rename(CellClass = Neuron_type1, SubType = Neuron_type3) %>%
  mutate(CellClass = recode(CellClass, "Excitatory" = "Exc", "Inhibitory" = "Inh"))

## compute density and fit parameters
sub_luo2017 <- do.call(rbind, map(1:nrow(info), function(i) wrapper2(info[i])))
fwrite(sub_luo2017, "data/interim/estim_emiBetaPrior_IBregression/emiBetaPrior_subtype_subsample_luo2017.csv")

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
sub_unmeth_smr0 <- sub_unmeth %>% 
  group_by(DataSource, CellClass, SubType, med_cov) %>% 
  summarise(weight = sum(cell_MF==0) / length(cell_MF), 
            mean = mean(cell_MF[cell_MF>0]), 
            var = var(cell_MF[cell_MF>0]))

# MOM estimates (plug in mean and var) vs. median coverage
sub_unmeth_smr0 %>% ggplot(aes(med_cov, weight)) + geom_point() + ylab("MOM estimate of nu") + xlab("Median coverage")
ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_paramNuMOM_vs_medCov.png", width = 5, height = 5)
sub_unmeth_smr0 %>% ggplot(aes(med_cov, mean)) + geom_point() + ylab("MOM estimate of mu") + xlab("Median coverage")
ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_paramMuMOM_vs_medCov.png", width = 5, height = 5)
sub_unmeth_smr0 %>% ggplot(aes(med_cov, ((mean*(1-mean)/var - 1)))) + geom_point()  + ylab("MOM estimate of sigma") + xlab("Median coverage")
ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_paramSigmaMOM_vs_medCov.png", width = 5, height = 5)
# (log-transform med_cov and suitable link functions on parameters)
sub_unmeth_smr0 %>% ggplot(aes(log(med_cov), weight)) + geom_point() + geom_smooth(method = "lm") + ylab("MOM estimate of nu") + xlab("Log(median coverage)")
ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_transformed_paramNuMOM_vs_medCov.png", width = 5, height = 5)
sub_unmeth_smr0 %>% ggplot(aes(log(med_cov), mean)) + geom_point() + geom_smooth(method = "lm") + ylab("MOM estimate of mu") + xlab("Log(median coverage)")
ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_transformed_paramMuMOM_vs_medCov.png", width = 5, height = 5)
sub_unmeth_smr0 %>% ggplot(aes(log(med_cov), ((mean*(1-mean)/var - 1))^(-1))) + geom_point() + geom_smooth(method = "lm") + ylab("Inverse MOM estimate of sigma") + xlab("Log(median coverage)")
ggsave("plots/estim_emiBetaPrior_IBregression/point_unmethPop_transformed_paramSigmaMOM_vs_medCov.png", width = 5, height = 5)

show.link(family = "BEZI")
mod_u <- gamlss(cell_MF ~ log(med_cov), data = sub_unmeth, 
              sigma.formula = ~ log(med_cov), nu.formula = ~ log(med_cov), 
              family=BEZI(nu.link = "identity", mu.link = "identity", sigma.link = "inverse"))

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
sub_meth_smr0 <- sub_meth %>% 
  group_by(DataSource, CellClass, SubType, med_cov) %>% 
  summarise(weight = sum(cell_MF==1) / length(cell_MF), 
            mean = mean(cell_MF[cell_MF<1]), 
            var = var(cell_MF[cell_MF<1]))

# MOM estimates (plug in mean and var) vs. median coverage
sub_meth_smr0 %>% ggplot(aes(med_cov, weight)) + geom_point() + ylab("MOM estimate of nu") + xlab("Median coverage")
ggsave("plots/estim_emiBetaPrior_IBregression/point_methPop_paramNuMOM_vs_medCov.png", width = 5, height = 5)
sub_meth_smr0 %>% ggplot(aes(med_cov, mean)) + geom_point() + ylab("MOM estimate of mu") + xlab("Median coverage")
ggsave("plots/estim_emiBetaPrior_IBregression/point_methPop_paramMuMOM_vs_medCov.png", width = 5, height = 5)
sub_meth_smr0 %>% ggplot(aes(med_cov, ((mean*(1-mean)/var - 1)))) + geom_point()  + ylab("MOM estimate of sigma") + xlab("Median coverage")
ggsave("plots/estim_emiBetaPrior_IBregression/point_methPop_paramSigmaMOM_vs_medCov.png", width = 5, height = 5)
# (log-transform med_cov and suitable link functions on parameters)
sub_meth_smr0 %>% ggplot(aes(log(med_cov), weight)) + geom_point() + geom_smooth(method = "lm") + ylab("MOM estimate of nu") + xlab("Log(median coverage)")
ggsave("plots/estim_emiBetaPrior_IBregression/point_methPop_transformed_paramNuMOM_vs_medCov.png", width = 5, height = 5)
sub_meth_smr0 %>% ggplot(aes(log(med_cov), mean)) + geom_point() + geom_smooth(method = "lm") + ylab("MOM estimate of mu") + xlab("Log(median coverage)")
ggsave("plots/estim_emiBetaPrior_IBregression/point_methPop_transformed_paramMuMOM_vs_medCov.png", width = 5, height = 5)
sub_meth_smr0 %>% ggplot(aes(log(med_cov), ((mean*(1-mean)/var - 1))^(-1))) + geom_point() + geom_smooth(method = "lm") + ylab("Inverse MOM estimate of sigma") + xlab("Log(median coverage)")
ggsave("plots/estim_emiBetaPrior_IBregression/point_methPop_transformed_paramSigmaMOM_vs_medCov.png", width = 5, height = 5)

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

