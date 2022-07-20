## Objective of this script: 
## To estimate emission beta-mixture priors from 18 cell subtypes from Liu2021 (QC'ed) data and 9 from (QC'ed) Luo2017 data. 
library(tidyverse)
setwd(here::here())
source("code/estim_emiBetaPrior2/helper_functions.R")

################################################################
# ==== Compute densities and fit parameters in Liu2021 data ====
################################################################

wrapper1 <- function(info, cores = 16){
  dir <- "data/processed/processed_liu2021/"
  file_list <- list.files(dir)
  file_name <- grep(pattern = paste0("sample", info$Sample,
                                     "_", info$GEO_accession, 
                                     "_subtype_", sub(" ", "_", info$SubType),
                                     ".*_qced$"),
                    x = file_list, value = T)
  print(file_name)
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))
  
  ### Changed here 
  dens <- calDensity(values(cells.se)$cell_MF, n_bins = 200)
  params_u <- calBestParamFromDensity(dens, population = "u", w1_space = seq(0, 1, 0.01), param2_space = seq(20,3000,10), mc.cores = cores)
  params_m <- calBestParamFromDensity(dens, population = "m", w1_space = seq(0, 1, 0.01), param2_space = seq(20,3000,10), mc.cores = cores)
  # dens <- calDensity(values(cells.se)$cell_MF, n_bins = 200)
  # params_u <- calBestParamFromDensity(dens, population = "u", param2_space = seq(1000,3000,10), mc.cores = cores)
  # params_m <- calBestParamFromDensity(dens, population = "m", param2_space = seq(1000,3000,10), mc.cores = cores)
  
  summary  <- data.table(Density = dens, 
                         Sample = info$Sample, 
                         GEO_accession = info$GEO_accession, 
                         CellClass = info$CellClass,
                         SubType = info$SubType, 
                         N_cell = ncol(cells.se),
                         w1_u = params_u$w1,
                         param1_u = params_u$param1,
                         param2_u = params_u$param2,
                         l1_diff_u = params_u$l1_diff[[1]],
                         w1_m = params_m$w1,
                         param1_m = params_m$param1,
                         param2_m = params_m$param2,
                         l1_diff_m = params_m$l1_diff[[1]]
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
fwrite(smr_liu2021, "data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_liu2021.csv")

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
  
  ### changed here
  dens <- calDensity(values(cells.se)$cell_MF, n_bins = 200)
  params_u <- calBestParamFromDensity(dens, population = "u", w1_space = seq(0, 1, 0.01), param2_space = seq(20,3000,10), mc.cores = cores)
  params_m <- calBestParamFromDensity(dens, population = "m", w1_space = seq(0, 1, 0.01), param2_space = seq(20,3000,10), mc.cores = cores)
  # dens <- calDensity(values(cells.se)$cell_MF, n_bins = 200)
  # params_u <- calBestParamFromDensity(dens, population = "u", param2_space = seq(1000,3000,10), mc.cores = cores)
  # params_m <- calBestParamFromDensity(dens, population = "m", param2_space = seq(1000,3000,10), mc.cores = cores)
  
  summary  <- data.table(Density = dens, 
                         CellClass = info$CellClass,
                         SubType = info$SubType, 
                         N_cell = ncol(cells.se),
                         w1_u = params_u$w1,
                         param1_u = params_u$param1,
                         param2_u = params_u$param2,
                         l1_diff_u = params_u$l1_diff[[1]],
                         w1_m = params_m$w1,
                         param1_m = params_m$param1,
                         param2_m = params_m$param2,
                         l1_diff_m = params_m$l1_diff[[1]]
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
fwrite(smr_luo2017, "data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_luo2017.csv")

################################################################
# ==== Fit linear model on parameters ~ N cell ====
################################################################
smr_liu2021 <- fread("data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_liu2021.csv") %>%
  select(-c(Sample, GEO_accession)) %>%
  mutate(DataSource = "Liu2021")
smr_luo2017 <- fread("data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_luo2017.csv") %>%
  mutate(DataSource = "Luo2017")

all(colnames(smr_liu2021)==colnames(smr_luo2017)) # check if the colnames are matched

smr_long <- rbind(smr_liu2021, smr_luo2017) # each 201-row is a subtype
smr_short <- smr_long %>%
  group_by(CellClass, SubType, N_cell, w1_u, param1_u, param2_u, l1_diff_u, w1_m, param1_m, param2_m, l1_diff_m) %>%
  summarise(DataSource = max(DataSource)) # each row is a subtype

fwrite(smr_long, "data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_with_subtypePDF.csv")
fwrite(smr_short, "data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_without_subtypePDF.csv")

smr_short <- fread("data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_without_subtypePDF.csv")
smr_long <- fread("data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_with_subtypePDF.csv")

## fitted parameters for unmeth subpop vs. N cells
lm_fit <- lm(w1_u ~ N_cell, data = smr_short)
plot(lm_fit, which = 4)
(lm <- summary(lm_fit))
smr_short %>%
  ggplot(aes(N_cell, w1_u)) +
  # geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  # annotate('text', x = 250, y = 0.99, label = paste("italic(R)^{2}==", round(lm$r.squared, 2)), parse = TRUE) +
  # annotate("text", x = 250, y = 0.95, label = paste0('p-value = ', signif(lm$coef[2,4], 2))) +
  ylab("weight") + xlab("N cells") + ylim(0, 1)
ggsave("plots/estim_emiBetaPrior2/point_fittedPrior_unmethPop_w1_vs_Ncell_27subtypes.png", width = 6, height = 5)

lm_fit <- lm(param1_u ~ N_cell, data = smr_short)
plot(lm_fit, which = 4)
(lm <- summary(lm_fit))
smr_short %>%
  ggplot(aes(N_cell, param1_u)) +
  # geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  # annotate("text", x = 250, y = 19, label = paste0('italic(R)^{2}==', round(lm$r.squared, 2)), parse = TRUE) +
  # annotate("text", x = 250, y = 18, label = paste0('p-value = ', signif(lm$coef[2,4], 2))) +
  ylab(expression(beta[1])) + xlab("N cells") + ylim(0,20)
ggsave("plots/estim_emiBetaPrior2/point_fittedPrior_unmethPop_param1_vs_Ncell_27subtypes.png", width = 6, height = 5)

lm_fit <- lm(param2_u ~ N_cell, data = smr_short)
plot(lm_fit, which = 4)
(lm <- summary(lm_fit))
smr_short %>%
  ggplot(aes(N_cell, param2_u)) +
  # geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) + 
  # annotate('text', x = 250, y = 2900, label = paste("italic(R)^{2}==", round(lm$r.squared, 2)), parse = TRUE) +
  # annotate("text", x = 250, y = 2780, label = paste0('p-value = ', signif(lm$coef[2,4], 2))) +
  ylab(expression(beta[2])) + xlab("N cells") + ylim(0,3000)
ggsave("plots/estim_emiBetaPrior2/point_fittedPrior_unmethPop_param2_vs_Ncell_27subtypes.png", width = 6, height = 5)

smr_short %>%
  ggplot(aes(N_cell, l1_diff_u)) +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  ylab("L1 difference") + xlab("N cells") + ylim(0, 1)
ggsave("plots/estim_emiBetaPrior2/point_fittedPrior_unmethPop_L1Diff_vs_Ncell_27subtypes.png", width = 6, height = 5)

## fitted parameters for meth subpop vs. N cells
lm_fit <- lm(w1_m ~ N_cell, data = smr_short)
plot(lm_fit, which = 4)
(lm <- summary(lm_fit))
smr_short %>%
  ggplot(aes(N_cell, w1_m)) +
  # geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) + 
  # annotate('text', x = 250, y = 0.99, label = paste("italic(R)^{2}==", round(lm$r.squared, 2)), parse = TRUE) +
  # annotate("text", x = 250, y = 0.95, label = paste0('p-value = ', signif(lm$coef[2,4], 2))) +
  ylab("weight") + xlab("N cells") + ylim(0, 1)
ggsave("plots/estim_emiBetaPrior2/point_fittedPrior_methPop_w1_vs_Ncell_27subtypes.png", width = 6, height = 5)

lm_fit <- lm(param1_m ~ N_cell, data = smr_short)
plot(lm_fit, which = 4)
(lm <- summary(lm_fit))
smr_short %>%
  ggplot(aes(N_cell, param1_m)) +
  # geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2)  + 
  # annotate("text", x = 250, y = 19, label = paste0('italic(R)^{2}==', round(lm$r.squared, 2)), parse = TRUE) +
  # annotate("text", x = 250, y = 18, label = paste0('p-value = ', signif(lm$coef[2,4], 2))) +
  ylab(expression(alpha[1])) + xlab("N cells") + ylim(0,20)
ggsave("plots/estim_emiBetaPrior2/point_fittedPrior_methPop_param1_vs_Ncell_27subtypes.png", width = 6, height = 5)

lm_fit <- lm(param2_m ~ N_cell, data = smr_short)
plot(lm_fit, which = 4)
(lm <- summary(lm_fit))
smr_short %>%
  ggplot(aes(N_cell, param2_m)) +
  # geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) + 
  # annotate('text', x = 250, y = 2900, label = paste("italic(R)^{2}==", round(lm$r.squared, 2)), parse = TRUE) +
  # annotate("text", x = 250, y = 2780, label = paste0('p-value = ', signif(lm$coef[2,4], 2))) +
  ylab(expression(alpha[2])) + xlab("N cells") + ylim(0,3000)
ggsave("plots/estim_emiBetaPrior2/point_fittedPrior_methPop_param2_vs_Ncell_27subtypes.png", width = 6, height = 5)

smr_short %>%
  ggplot(aes(N_cell, l1_diff_m)) +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) + 
  ylab("L1 difference") + xlab("N cells") + ylim(0, 1)
ggsave("plots/estim_emiBetaPrior2/point_fittedPrior_methPop_L1Diff_vs_Ncell_27subtypes.png", width = 6, height = 5)


########################
# ==== FINAL MODELS ====
########################

## unmeth subpop:
smr_short <- fread("data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_without_subtypePDF.csv")
# lm(w1_u ~ N_cell, data = smr_short) # w
# lm(param1_u ~ N_cell, data = smr_short) # beta_1 (stops changing once w1_u hits 1)
# lm(param2_u ~ 1, data = smr_short) # beta_2

## meth subpop:
smr_short <- fread("data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_without_subtypePDF.csv")
# lm(w1_m ~ N_cell, data = smr_short) # w
# lm(param1_m ~ N_cell, data = smr_short) # beta_1 (stops changing once w1_u hits 1)
# lm(param2_m ~ 1, data = smr_short) # beta_2

################################################################
# ==== Predict parameters with fitted linear model ====
################################################################

smr_short <- fread("data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_without_subtypePDF.csv")
smr_long <- fread("data/interim/estim_emiBetaPrior2/emiBetaPrior_subtype_summary_with_subtypePDF.csv")

### import function with final prediction model incorporated
.priorParams <- function(x) {0}
insertSource("code/package_functions/helper_functions.R", functions=".priorParams") 

smr_long_pred <- smr_short %>%
  mutate(w1_u_pred = map_dbl(N_cell, ~.priorParams(.x, type = "u")[1]),
         param1_u_pred = map_dbl(N_cell, ~.priorParams(.x, type = "u")[2]),
         param2_u_pred = map_dbl(N_cell, ~.priorParams(.x, type = "u")[3]),
         w1_m_pred = map_dbl(N_cell, ~.priorParams(.x, type = "m")[1]),
         param1_m_pred = map_dbl(N_cell, ~.priorParams(.x, type = "m")[2]),
         param2_m_pred = map_dbl(N_cell, ~.priorParams(.x, type = "m")[3])) %>% 
  right_join(smr_long)

## compute L1 difference between true data distribution and predicted distribution from linear model
#### in unmethylated subpopulation
smr_long_pred_u <- smr_long_pred %>%
  select(-c(param1_m, param1_m_pred, param2_m, param2_m_pred, w1_m, w1_m_pred, l1_diff_m)) %>%
  group_by(SubType, N_cell) %>%
  mutate(Density.y = replace(Density.y, Density.x >= 0.5, 0)) %>%
  mutate(Density.y = Density.y / sum(Density.y)) %>% # normalize
  mutate(Density.y_pred = w1_u*dbeta(Density.x, 1, param1_u) + (1-w1_u) * dbeta(Density.x, 1, param2_u)) %>% # predicted density
  mutate(Density.y_pred = Density.y_pred / sum(Density.y_pred)) %>% # normalize
  mutate(l1_diff_u_pred = sum(abs(Density.y_pred*1e5-Density.y*1e5))/1e5) # compute l1 difference

smr_short_pred_u <- smr_long_pred_u %>%
  select(-c(Density.x, Density.y, Density.y_pred)) %>%
  dplyr::slice(seq(1, nrow(smr_long_pred_u), 201))

smr_short_pred_u %>%
  ggplot(aes(l1_diff_u, l1_diff_u_pred)) +
  geom_point(aes(color = DataSource, shape = CellClass, size = N_cell)) + 
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  xlab("L1 difference (original)") + ylab("L1 difference (predicted)") + 
  xlim(0,1.1) + ylim(0, 1.1)
ggsave("plots/estim_emiBetaPrior2/point_unmethPop_predicted_vs_original_L1Diff_27subtypes.png", width = 6, height = 5)

#### in methylated subpopulation
smr_long_pred_m <- smr_long_pred %>%
  select(-c(param1_u, param1_u_pred, param2_u, param2_u_pred, w1_u, w1_u_pred, l1_diff_u)) %>%
  group_by(SubType, N_cell) %>%
  mutate(Density.y = replace(Density.y, Density.x <= 0.5, 0)) %>%
  mutate(Density.y = Density.y / sum(Density.y)) %>% # normalize
  mutate(Density.y_pred = w1_m*dbeta(Density.x, param1_m, 1) + (1-w1_m) * dbeta(Density.x, param2_m, 1)) %>% # predicted density
  mutate(Density.y_pred = Density.y_pred / sum(Density.y_pred)) %>% # normalize
  mutate(l1_diff_m_pred = sum(abs(Density.y_pred*1e5-Density.y*1e5))/1e5) # compute l1 difference

smr_short_pred_m <- smr_long_pred_m %>%
  select(-c(Density.x, Density.y, Density.y_pred)) %>%
  dplyr::slice(seq(1, nrow(smr_long_pred_m), 201))

smr_short_pred_m %>%
  ggplot(aes(l1_diff_m, l1_diff_m_pred)) +
  geom_point(aes(color = DataSource, shape = CellClass, size = N_cell)) + 
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  xlab("L1 difference (original)") + ylab("L1 difference (predicted)") + 
  xlim(0,1.1) + ylim(0, 1.1)
ggsave("plots/estim_emiBetaPrior2/point_methPop_predicted_vs_original_L1Diff_27subtypes.png", width = 6, height = 5)


