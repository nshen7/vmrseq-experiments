## Objective of this script: 
## To estimate emission beta-mixture priors from 18 cell subtypes from Liu2021 (QC'ed) data and 9 from (QC'ed) Luo2017 data. 
library(tidyverse)
setwd(here::here())


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
  
  mf_u <- values(cells.se)$cell_MF[values(cells.se)$cell_MF < 0.5]
  w_u <- sum(mf_u == 0) / length(mf_u)
  mean_u <- mean(mf_u[mf_u>0]); var_u <- var(mf_u[mf_u>0])

  mf_m <- values(cells.se)$cell_MF[values(cells.se)$cell_MF > 0.5]
  w_m <- sum(mf_m == 1) / length(mf_m)
  mean_m <- mean(mf_m[mf_m>0]); var_m <- var(mf_m[mf_m>0])
  
  summary  <- data.table(Sample = info_row$Sample, 
                         GEO_accession = info_row$GEO_accession, 
                         CellClass = info_row$CellClass,
                         SubType = info_row$SubType, 
                         N_cell = ncol(cells.se),
                         w_u = w_u,
                         mean_u = mean_u,
                         var_u = var_u,
                         alpha_u = (mean_u*(1-mean_u)/var_u - 1) * mean_u,
                         beta_u = (mean_u*(1-mean_u)/var_u - 1) * (1-mean_u),
                         w_m = w_m,
                         mean_m = mean_m,
                         var_m = var_m,
                         alpha_m = (mean_m*(1-mean_m)/var_m - 1) * mean_m,
                         beta_m = (mean_m*(1-mean_m)/var_m - 1) * (1-mean_m)
  )
  return(summary)
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
smr_liu2021 <- do.call(rbind, map(1:nrow(info), function(i) wrapper1(info[i])))
fwrite(smr_liu2021, "data/interim/estim_emiBetaPrior_inflated/emiBetaPrior_subtype_summary_liu2021.csv")


################################################################
# ==== Compute densities and fit parameters in Luo2017 data ====
################################################################

wrapper2 <- function(info, cores = 16){
  dir <- "data/processed/processed_luo2017_mice/"
  file_list <- list.files(dir)
  file_name <- grep(pattern = paste0("subtype_", sub("/", "", info$SubType), ".*_qced$"),
                    x = file_list, value = T)
  print(file_name)
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))
  
  mf_u <- values(cells.se)$cell_MF[values(cells.se)$cell_MF < 0.5]
  w_u <- sum(mf_u == 0) / length(mf_u)
  mean_u <- mean(mf_u[mf_u>0]); var_u <- var(mf_u[mf_u>0])
  
  mf_m <- values(cells.se)$cell_MF[values(cells.se)$cell_MF > 0.5]
  w_m <- sum(mf_m == 1) / length(mf_m)
  mean_m <- mean(mf_m[mf_m>0]); var_m <- var(mf_m[mf_m>0])

  summary  <- data.table(CellClass = info$CellClass,
                         SubType = info$SubType, 
                         N_cell = ncol(cells.se),
                         w_u = w_u,
                         mean_u = mean_u,
                         var_u = var_u,
                         alpha_u = (mean_u*(1-mean_u)/var_u - 1) * mean_u,
                         beta_u = (mean_u*(1-mean_u)/var_u - 1) * (1-mean_u),
                         w_m = w_m,
                         mean_m = mean_m,
                         var_m = var_m,
                         alpha_m = (mean_m*(1-mean_m)/var_m - 1) * mean_m,
                         beta_m = (mean_m*(1-mean_m)/var_m - 1) * (1-mean_m)
  )
  
  return(summary)
}

## information of all subtypes used for training
metadata <- fread("../../DXM_extend_chr1/data/metadata/sample_info_processed.csv")
info <- metadata[, .(.N), by = .(Neuron_type1, Neuron_type3)] %>% 
  arrange(desc(N)) %>%
  filter(N >= 100) %>%
  dplyr::rename(CellClass = Neuron_type1, SubType = Neuron_type3) %>%
  mutate(CellClass = recode(CellClass, "Excitatory" = "Exc", "Inhibitory" = "Inh"))

## compute density and fit parameters
smr_luo2017 <- do.call(rbind, map(1:nrow(info), function(i) wrapper2(info[i])))
fwrite(smr_luo2017, "data/interim/estim_emiBetaPrior_inflated/emiBetaPrior_subtype_summary_luo2017.csv")

##################################
# ==== Write out summary data ====
##################################

smr_liu2021 <- fread("data/interim/estim_emiBetaPrior_inflated/emiBetaPrior_subtype_summary_liu2021.csv") %>%
  dplyr::select(-c(Sample, GEO_accession)) %>%
  mutate(DataSource = "Liu2021")
smr_luo2017 <- fread("data/interim/estim_emiBetaPrior_inflated/emiBetaPrior_subtype_summary_luo2017.csv") %>%
  mutate(DataSource = "Luo2017")

all(colnames(smr_liu2021)==colnames(smr_luo2017)) # check if the colnames are matched

smr <- rbind(smr_liu2021, smr_luo2017)

fwrite(smr, "data/interim/estim_emiBetaPrior_inflated/emiBetaPrior_subtype_summary.csv")


###################################################
# ==== Fit linear model on parameters ~ N cell ====
###################################################

smr <- fread("data/interim/estim_emiBetaPrior_inflated/emiBetaPrior_subtype_summary.csv")

#### Add in 2 sample-combined subtypes
# dirs <- c("data/processed/processed_liu2021/summary_subtype_CA1_Chrm3_1359cells_combined_qced.rds",
#           "data/processed/processed_liu2021/summary_subtype_DG_dg-all_1462cells_combined_qced.rds")
# for (i in 1:2) {
#   df <- readRDS(dirs[i])
#   mf_u <- df$cell_MF[df$cell_MF < 0.5]
#   w_u <- sum(mf_u == 0) / length(mf_u)
#   mean_u <- mean(mf_u[mf_u>0]); var_u <- var(mf_u[mf_u>0])
#   
#   mf_m <- df$cell_MF[df$cell_MF > 0.5]
#   w_m <- sum(mf_m == 1) / length(mf_m)
#   mean_m <- mean(mf_m[mf_m>0]); var_m <- var(mf_m[mf_m>0])
#   
#   if(i==1) {
#     CellClass = "Exc"; SubType = "CA1 Chrm3"; N_cell = 1359
#   } else {
#     CellClass = "Exc"; SubType = "DG dg-all"; N_cell = 1462
#   }
#   summary  <- data.table(CellClass = CellClass,
#                          SubType = SubType, 
#                          N_cell = N_cell,
#                          w_u = w_u,
#                          mean_u = mean_u,
#                          var_u = var_u,
#                          alpha_u = (mean_u*(1-mean_u)/var_u - 1) * mean_u,
#                          beta_u = (mean_u*(1-mean_u)/var_u - 1) * (1-mean_u),
#                          w_m = w_m,
#                          mean_m = mean_m,
#                          var_m = var_m,
#                          alpha_m = (mean_m*(1-mean_m)/var_m - 1) * mean_m,
#                          beta_m = (mean_m*(1-mean_m)/var_m - 1) * (1-mean_m),
#                          DataSource = "Combined_liu"
#   )
#   summary
#   smr <- rbind(smr, summary)
# }

## fitted parameters for unmeth subpop vs. N cells
(loglm <- summary(lm(w_u ~ log(N_cell), data = smr)))
(lm <- summary(lm(w_u ~ N_cell, data = smr)))
smr %>%
  ggplot(aes(N_cell, w_u)) +
  geom_smooth(method = "lm", formula = y~log(x), se = F, size = 0.7, color = "skyblue3") +
  geom_smooth(method = "lm", se = F, linetype = "dashed", size = 0.7, color = "orange") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  annotate("text", x = 350, y = 0.02, label = paste0('Log-linear: R^2 = ', round(loglm$r.squared, 3)), color = "skyblue3") +
  annotate('text', x = 350, y = 0.06, fontface = 2, label = paste0("Linear: R^2 = ", round(lm$r.squared, 3)), color = "orange") +
  ylab("Prob mass at 0") + xlab("N cells") + ylim(0,1)
ggsave(paste0("plots/estim_emiBetaPrior_inflated/point_fittedPrior_unmethPop_w_vs_Ncell_",nrow(smr),"subtypes.png"), width = 6, height = 5)

(loglm <- summary(lm(alpha_u ~ log(N_cell), data = smr)))
(lm <- summary(lm(alpha_u ~ N_cell, data = smr)))
smr %>%
  ggplot(aes(N_cell, alpha_u)) +
  geom_smooth(method = "lm", formula = y~log(x), se = F, size = 0.7, color = "skyblue3") +
  geom_smooth(method = "lm", se = F, linetype = "dashed", size = 0.7, color = "orange") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  annotate("text", x = 350, y = 0.15, fontface = 2, label = paste0('Log-linear model: R^2 = ', round(loglm$r.squared, 3)), color = "skyblue3") +
  annotate('text', x = 350, y = 0.3, label = paste0("Linear model: R^2 = ", round(lm$r.squared, 3)), color = "orange") +
  ylab("Alpha") + xlab("N cells") + ylim(0, 3)
ggsave(paste0("plots/estim_emiBetaPrior_inflated/point_fittedPrior_unmethPop_alpha_vs_Ncell_",nrow(smr),"subtypes.png"), width = 6, height = 5)

(loglm <- summary(lm(beta_u ~ log(N_cell), data = smr)))
(lm <- summary(lm(beta_u ~ N_cell, data = smr)))
smr %>%
  ggplot(aes(N_cell, beta_u)) +
  geom_smooth(method = "lm", formula = y~log(x), se = F, size = 0.7, color = "skyblue3") +
  geom_smooth(method = "lm", se = F, linetype = "dashed", size = 0.7, color = "orange") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  annotate("text", x = 350, y = 3, fontface = 2, label = paste0('Log-linear model: R^2 = ', round(loglm$r.squared, 3)), color = "skyblue3") +
  annotate('text', x = 350, y = 3.3, label = paste0("Linear model: R^2 = ", round(lm$r.squared, 3)), color = "orange") +
  ylab("Beta") + xlab("N cells") + ylim(3, 9)
ggsave(paste0("plots/estim_emiBetaPrior_inflated/point_fittedPrior_unmethPop_beta_vs_Ncell_",nrow(smr),"subtypes.png"), width = 6, height = 5)

## fitted parameters for meth subpop vs. N cells
(loglm <- summary(lm(w_m ~ log(N_cell), data = smr)))
(lm <- summary(lm(w_m ~ N_cell, data = smr)))
smr %>%
  ggplot(aes(N_cell, w_m)) +
  geom_smooth(method = "lm", formula = y~log(x), se = F, size = 0.7, color = "skyblue3") +
  geom_smooth(method = "lm", se = F, linetype = "dashed", size = 0.7, color = "orange") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  annotate("text", x = 300, y = 0.02, fontface = 2, label = paste0('Log-linear: R^2 = ', round(loglm$r.squared, 3)), color = "skyblue3") +
  annotate('text', x = 300, y = 0.06, label = paste0("Linear: R^2 = ", round(lm$r.squared, 3)), color = "orange") +
  ylab("Prob mass at 1") + xlab("N cells") + ylim(0, 1)
ggsave(paste0("plots/estim_emiBetaPrior_inflated/point_fittedPrior_methPop_w_vs_Ncell_",nrow(smr),"subtypes.png"), width = 6, height = 5)

(loglm <- summary(lm(alpha_m ~ log(N_cell), data = smr)))
(lm <- summary(lm(alpha_m ~ N_cell, data = smr)))
smr %>%
  ggplot(aes(N_cell, alpha_m)) +
  geom_smooth(method = "lm", formula = y~log(x), se = F, size = 0.7, color = "skyblue3") +
  geom_smooth(method = "lm", se = F, linetype = "dashed", size = 0.7, color = "orange") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  annotate("text", x = 350, y = 7, fontface = 2, label = paste0('Log-linear model: R^2 = ', round(loglm$r.squared, 3)), color = "skyblue3") +
  annotate('text', x = 350, y = 7.2, label = paste0("Linear model: R^2 = ", round(lm$r.squared, 3)), color = "orange") +
  ylab("Alpha") + xlab("N cells") + ylim(3.8, 7.2)
ggsave(paste0("plots/estim_emiBetaPrior_inflated/point_fittedPrior_methPop_alpha_vs_Ncell_",nrow(smr),"subtypes.png"), width = 6, height = 5)

(loglm <- summary(lm(beta_m ~ log(N_cell), data = smr)))
(lm <- summary(lm(beta_m ~ N_cell, data = smr)))
smr %>%
  ggplot(aes(N_cell, beta_m)) +
  geom_smooth(method = "lm", formula = y~log(x), se = F, size = 0.7, color = "skyblue3") +
  geom_smooth(method = "lm", se = F, linetype = "dashed", size = 0.7, color = "orange") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  annotate("text", x = 350, y = 0.7, fontface = 2, label = paste0('Log-linear model: R^2 = ', round(loglm$r.squared, 3)), color = "skyblue3") +
  annotate('text', x = 350, y = 0.73, label = paste0("Linear model: R^2 = ", round(lm$r.squared, 3)), color = "orange") +
  ylab("Beta") + xlab("N cells") + ylim(0.25,0.75)
ggsave(paste0("plots/estim_emiBetaPrior_inflated/point_fittedPrior_methPop_beta_vs_Ncell_",nrow(smr),"subtypes.png"), width = 6, height = 5)

smr %>%
  ggplot(aes(N_cell, alpha_u/(alpha_u+beta_u))) +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) 
smr %>%
  ggplot(aes(N_cell, alpha_m/(alpha_m+beta_m))) +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) 
smr %>%
  ggplot(aes(N_cell, alpha_u*beta_u/(alpha_u+beta_u)^2/(alpha_u+beta_u+1))) +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) 
smr %>%
  ggplot(aes(N_cell, alpha_m/(alpha_m+beta_m))) +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2)


########################
# ==== FINAL MODELS ====
########################

## unmeth subpop:
smr <- fread("data/interim/estim_emiBetaPrior_inflated/emiBetaPrior_subtype_summary.csv")
# lm(w1_u ~ N_cell, data = smr) # w
# lm(alpha_u ~ N_cell, data = smr) # beta_1 (stops changing once w1_u hits 1)
# lm(beta_u ~ 1, data = smr) # beta_2

## meth subpop:
smr <- fread("data/interim/estim_emiBetaPrior_inflated/emiBetaPrior_subtype_summary.csv")
# lm(w1_m ~ N_cell, data = smr) # w
# lm(alpha_m ~ N_cell, data = smr) # beta_1 (stops changing once w1_u hits 1)
# lm(beta_m ~ 1, data = smr) # beta_2

