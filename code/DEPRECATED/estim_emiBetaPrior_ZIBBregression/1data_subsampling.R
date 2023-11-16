suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(data.table))

setwd(here::here())

## number of sites sampled from each subtype
sample_n <- 50000
set.seed(2022)

################################################################
# ==== Subsample subtypes in Liu2021 mice data ====
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
  
  # sub_dat_glm <- sub_dat %>%
  #   mutate(y_u = as.matrix(data.frame(cell_meth, cell_cov-cell_meth)),
  #          y_m = as.matrix(data.frame(cell_cov-cell_meth, cell_meth))) 
  # mod_u <- gamlss(y_u ~ 1, data = sub_dat_glm %>% filter(cell_MF < 0.5), family = ZIBB, 
  #                 trace = F, n.cyc = 100)
  # mod_m <- gamlss(y_m ~ 1, data = sub_dat_glm %>% filter(cell_MF > 0.5), family = ZIBB, 
  #                 trace = F, n.cyc = 100)
  # 
  # sub_dat$mle_nu_u <- fitted(mod_u, "nu")[1]
  # sub_dat$mle_mu_u <- fitted(mod_u, "mu")[1]
  # sub_dat$mle_sigma_u <- unname(fitted(mod_u, "sigma")[1])
  # sub_dat$mle_nu_m <- fitted(mod_m, "nu")[1]
  # sub_dat$mle_mu_m <- fitted(mod_m, "mu")[1]
  # sub_dat$mle_sigma_m <- unname(fitted(mod_m, "sigma")[1])
  
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
fwrite(sub_liu2021, "data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_liu2021.csv")


################################################################
# ==== Subsample subtypes in Luo2017 mice data ====
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
    DataSource = "Luo2017m",
    CellClass = info_row$CellClass,
    SubType = info_row$SubType,
    N_cell = ncol(cells.se),
    med_cov = median(values(cells.se)$cell_cov),
    cell_cov = values(cells.se)$cell_cov[index],
    cell_meth = values(cells.se)$cell_meth[index],
    cell_MF = values(cells.se)$cell_MF[index]
  )
  
  # sub_dat_glm <- sub_dat %>%
  #   mutate(y_u = as.matrix(data.frame(cell_meth, cell_cov-cell_meth)),
  #          y_m = as.matrix(data.frame(cell_cov-cell_meth, cell_meth))) 
  # mod_u <- gamlss(y_u ~ 1, data = sub_dat_glm %>% filter(cell_MF < 0.5), family = ZIBB, 
  #                 trace = F, n.cyc = 100)
  # mod_m <- gamlss(y_m ~ 1, data = sub_dat_glm %>% filter(cell_MF > 0.5), family = ZIBB, 
  #                 trace = F, n.cyc = 100)
  # 
  # sub_dat$mle_nu_u <- fitted(mod_u, "nu")[1]
  # sub_dat$mle_mu_u <- fitted(mod_u, "mu")[1]
  # sub_dat$mle_sigma_u <- unname(fitted(mod_u, "sigma")[1])
  # sub_dat$mle_nu_m <- fitted(mod_m, "nu")[1]
  # sub_dat$mle_mu_m <- fitted(mod_m, "mu")[1]
  # sub_dat$mle_sigma_m <- unname(fitted(mod_m, "sigma")[1])
  
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
sub_luo2017m <- do.call(rbind, map(1:nrow(info), function(i) wrapper2(info[i])))
fwrite(sub_luo2017m, "data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_luo2017.csv")


################################################################
# ==== Subsample subtypes in Luo2017 human data ====
################################################################

# wrapper3 <- function(info_row, cores = 16){
#   dir <- "data/processed/processed_luo2017_human/"
#   file_list <- list.files(dir)
#   file_name <- grep(pattern = paste0("subtype_", sub("/", "", info_row$SubType), ".*_qced$"),
#                     x = file_list, value = T)
#   print(file_name)
#   cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))
#   
#   index <- sample(1:nrow(cells.se), sample_n)
#   sub_dat <- data.frame(
#     DataSource = "Luo2017h",
#     CellClass = "NA",
#     SubType = info_row$SubType,
#     N_cell = ncol(cells.se),
#     med_cov = median(values(cells.se)$cell_cov),
#     cell_cov = values(cells.se)$cell_cov[index],
#     cell_meth = values(cells.se)$cell_meth[index],
#     cell_MF = values(cells.se)$cell_MF[index]
#   )
#   
#   sub_dat_glm <- sub_dat %>%
#     mutate(y_u = as.matrix(data.frame(cell_meth, cell_cov-cell_meth)),
#            y_m = as.matrix(data.frame(cell_cov-cell_meth, cell_meth))) 
#   mod_u <- gamlss(y_u ~ 1, data = sub_dat_glm %>% filter(cell_MF < 0.5), family = ZIBB, 
#                   trace = F, n.cyc = 100)
#   mod_m <- gamlss(y_m ~ 1, data = sub_dat_glm %>% filter(cell_MF > 0.5), family = ZIBB, 
#                   trace = F, n.cyc = 100)
#   
#   sub_dat$mle_nu_u <- fitted(mod_u, "nu")[1]
#   sub_dat$mle_mu_u <- fitted(mod_u, "mu")[1]
#   sub_dat$mle_sigma_u <- unname(fitted(mod_u, "sigma")[1])
#   sub_dat$mle_nu_m <- fitted(mod_m, "nu")[1]
#   sub_dat$mle_mu_m <- fitted(mod_m, "mu")[1]
#   sub_dat$mle_sigma_m <- unname(fitted(mod_m, "sigma")[1])
#   
#   return(sub_dat)
# } 
# 
# ## information of all subtypes used for training
# metadata <- fread("data/metadata/metadata_luo2017/NIHMS893063-supplement-Table_S2_csv.csv", skip = 1)[, .(Sample, `Neuron type`)]
# info <- metadata[, .(.N), by = .(`Neuron type`)] %>%
#   arrange(desc(N)) %>%
#   filter(N >= 140) %>%
#   dplyr::rename(SubType = `Neuron type`) 
# 
# ## compute density and fit parameters
# sub_luo2017h <- do.call(rbind, map(1:nrow(info), function(i) wrapper3(info[i])))

#############################################
#### Processed into meth and unmeth data ####
#############################################

sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_liu2021.csv")
sub_luo2017 <- fread("data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_luo2017.csv")
sub_unmeth <- rbind(sub_luo2017, sub_liu2021) %>% 
  as_tibble() %>%
  dplyr::filter(cell_MF < 0.5) %>%
  mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth)))
saveRDS(sub_unmeth, "data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_unmethClust.rds")
sub_meth <- rbind(sub_luo2017, sub_liu2021) %>% 
  as_tibble() %>%
  dplyr::filter(cell_MF > 0.5) %>%
  mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth)))
saveRDS(sub_meth, "data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_methClust.rds")

#### Optional: Add in 2 sample-combined subtypes
# dirs <- c("data/processed/processed_liu2021/summary_subtype_CA1_Chrm3_1359cells_combined_qced.rds",
#           "data/processed/processed_liu2021/summary_subtype_DG_dg-all_1462cells_combined_qced.rds")
# for (i in 1:2) {
#   df <- readRDS(dirs[i])
#   set.seed(2022); index <- sample(1:nrow(df), sample_n)
#   sub_dat <- data.frame(
#     DataSource = "Liu2021",
#     CellClass = "Exc",
#     SubType = ifelse(i == 1, "CA1 Chrm3", "DG dg-all"),
#     N_cell = ifelse(i == 1, 1359, 1462),
#     med_cov = median(df$cell_cov),
#     cell_cov = df$cell_cov[index],
#     cell_meth = df$cell_meth[index],
#     cell_MF = df$cell_MF[index]
#   )
#   
#   sub_dat_glm <- sub_dat %>%
#     mutate(y_u = as.matrix(data.frame(cell_meth, cell_cov-cell_meth)),
#            y_m = as.matrix(data.frame(cell_cov-cell_meth, cell_meth)))
#   mod_u <- gamlss(y_u ~ 1, data = sub_dat_glm %>% filter(cell_MF < 0.5), family = ZIBB,
#                   trace = F, n.cyc = 100)
#   mod_m <- gamlss(y_m ~ 1, data = sub_dat_glm %>% filter(cell_MF > 0.5), family = ZIBB,
#                   trace = F, n.cyc = 100)
# 
#   sub_dat$mle_nu_u <- fitted(mod_u, "nu")[1]
#   sub_dat$mle_mu_u <- fitted(mod_u, "mu")[1]
#   sub_dat$mle_sigma_u <- unname(fitted(mod_u, "sigma")[1])
#   sub_dat$mle_nu_m <- fitted(mod_m, "nu")[1]
#   sub_dat$mle_mu_m <- fitted(mod_m, "mu")[1]
#   sub_dat$mle_sigma_m <- unname(fitted(mod_m, "sigma")[1])
#   
#   sub_liu2021 <- rbind(sub_liu2021, sub_dat)
# }

# sub_unmeth <- rbind(sub_luo2017h, sub_luo2017m, sub_liu2021) %>%
#   as_tibble() %>%
#   dplyr::filter(cell_MF < 0.5) %>%
#   dplyr::select(-c(mle_nu_m, mle_mu_m, mle_sigma_m)) %>%
#   mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth)))
# sub_meth <- rbind(sub_luo2017h, sub_luo2017m, sub_liu2021) %>%
#   as_tibble() %>%
#   dplyr::filter(cell_MF > 0.5) %>%
#   dplyr::select(-c(mle_nu_u, mle_mu_u, mle_sigma_u)) %>%
#   mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth)))

# smr_unmeth <- sub_unmeth %>%
#   group_by(SubType, DataSource, CellClass, med_cov) %>%
#   summarise(mle_nu_u = max(mle_nu_u), mle_mu_u = max(mle_mu_u), mle_sigma_u = max(mle_sigma_u))
# smr_unmeth %>%
#   ggplot(aes(med_cov, mle_nu_u, color = DataSource, shape = CellClass)) +
#   geom_point()
# smr_unmeth %>%
#   ggplot(aes(med_cov, mle_mu_u, color = DataSource, shape = CellClass)) +
#   geom_point()
# smr_unmeth %>%
#   ggplot(aes(med_cov, mle_sigma_u, color = DataSource, shape = CellClass)) +
#   geom_point()
# 
# smr_meth <- sub_meth %>% 
#   group_by(SubType, DataSource, CellClass, med_cov) %>%
#   summarise(mle_nu_m = max(mle_nu_m), mle_mu_m = max(mle_mu_m), mle_sigma_m = max(mle_sigma_m)) 
# smr_meth %>%
#   ggplot(aes(med_cov, mle_nu_m, color = DataSource, shape = CellClass)) +
#   geom_point()
# smr_meth %>%
#   ggplot(aes(med_cov, mle_mu_m, color = DataSource, shape = CellClass)) +
#   geom_point()
# smr_meth %>%
#   ggplot(aes(med_cov, mle_sigma_m, color = DataSource, shape = CellClass)) +
#   geom_point()

