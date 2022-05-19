library(SummarizedExperiment)
library(GenomicRanges)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
library(bsseq)
library(tidyverse)
setAutoBlockSize(size=1e9)
setwd(here::here())
source("code/estim_emiBetaPrior/helper_functions.R")

n_bins <- 200

# ==== read in metadata ====
metadata <- fread("../../DXM_extend_chr1/data/metadata/sample_info_processed.csv")
cell_counts <- metadata[, .(.N), by = .(Neuron_type1, Neuron_type3)] %>% arrange(desc(N)) 

dt_info <- cell_counts[Neuron_type3 %in% c("mL6-2", "mL2/3", "mL4", "mDL-2", "mPv", "mSst-12")]
dt_info
#    Neuron_type1 Neuron_type3   N
# 1:   Excitatory        mL6-2 686
# 2:   Excitatory        mL2/3 649
# 3:   Excitatory          mL4 370
# 4:   Excitatory        mDL-2 272
# 5:   Inhibitory          mPv 136
# 6:   Inhibitory      mSst-12 123

dt_dens <- data.table()
for (i in 1:nrow(dt_info)) {  
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_mice/subtype_", sub("/", "", dt_info$Neuron_type3[i]), 
                                                        "_", dt_info$N[i], "cells_qced"))
  dens <- calDensity(values(cells.se)$cell_MF, n_bins = n_bins)
  params_u <- calBestParamFromDensity(dens, population = "u", mc.cores = 16)
  params_m <- calBestParamFromDensity(dens, population = "m", mc.cores = 16)

    dt_dens <- rbind(dt_dens, data.table(Density = dens,
                                         CellClass = dt_info$Neuron_type1[i],
                                         SubType = dt_info$Neuron_type3[i], 
                                         N_cell = dt_info$N[i],
                                         MedCovPct = metadata[Neuron_type3==dt_info$Neuron_type3[i], median(`Coverage(%)`)],
                                         w1_u = params_u$w1,
                                         param1_u = params_u$param1,
                                         param2_u = params_u$param2,
                                         l1_diff_u = params_u$l1_diff[[1]],
                                         w1_m = params_m$w1,
                                         param1_m = params_m$param1,
                                         param2_m = params_m$param2,
                                         l1_diff_m = params_m$l1_diff[[1]]
                                         )
                     )
  print(i)
}

dt_dens <- dt_dens %>%
  group_by(SubType) %>%
  mutate(CumDens = cumsum(Density.y))

write.table(dt_dens, "data/interim/estim_emiBetaPrior/eda/estimated_density_Luo2017_mice_6subtypes.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")

dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Luo2017_mice_6subtypes.txt") 
dt_dens %>% ggplot(aes(x = Density.x, y = CumDens)) + 
  geom_path(aes(linetype = CellClass, group = SubType, color = SubType)) +
  xlab("Value") +
  ylab("CDF of across-cell MF from subtypes in Luo2017 mice data") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_Luo2017_mice_6subtypes_coloredBySubtype.png", width = 7, height = 5)

dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Luo2017_mice_6subtypes.txt") 
dt_dens %>% ggplot(aes(x = Density.x, y = CumDens)) + 
  geom_path(aes(linetype = CellClass, group = SubType, color = MedCovPct)) +
  viridis::scale_color_viridis(name = "Median coverage (%)", option = "D") +
  xlab("Value") +
  ylab("CDF of across-cell MF from subtypes in Luo2017 mice data") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_Luo2017_mice_6subtypes_coloredByCovPct.png", width = 7, height = 5)

dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Luo2017_mice_6subtypes.txt") 
dt_dens %>% ggplot(aes(x = Density.x, y = CumDens)) + 
  geom_path(aes(linetype = CellClass, group = SubType, color = N_cell)) +
  viridis::scale_color_viridis(name = "N cells", option = "D") +
  xlab("Value") +
  ylab("CDF of across-cell MF from subtypes in Luo2017 mice data") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_Luo2017_mice_6subtypes_coloredByNcell.png", width = 7, height = 5)

dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Luo2017_mice_6subtypes.txt") 
dt_dens %>%
  filter(SubType=="mSst-12") %>%
  ggplot(aes(Density.x, Density.y)) +
  geom_bar(stat = "identity")



# # ==== estimate beta-mixture prior parameters ====
# best_params_u <- calBestParam(values(cells.se)$cell_MF, population = "u", mc.cores = 16)
# # Fitted mixture distribution density for unmethylated population: 0.84 * dbeta(x, 1, 4) + 0.16 * dbeta(x, 1, 980) 
# # L1 difference between discrete empirical distribution and fitted mixture:0.214801532557046
# best_params_m <- calBestParam(values(cells.se)$cell_MF, population = "m", mc.cores = 16)
# # Fitted mixture distribution density for methylated population: 0.96 * dbeta(x, 14, 1) + 0.04 * dbeta(x, 910, 1) 
# # L1 difference between discrete empirical distribution and fitted mixture:0.350341513900654
# 
# x <- seq(from = 0, to = 1, by = 1/n_bins)
# u <- 0.84 * dbeta(x, 1, 4) + 0.16 * dbeta(x, 1, 980) 
# m <- 0.96 * dbeta(x, 14, 1) + 0.04 * dbeta(x, 910, 1) 
# u[1]; m[n_bins+1] # density on 0 and 1
# 
# hist(values(cells.se)$cell_MF[values(cells.se)$cell_MF<0.4], breaks = 0.4*n_bins, prob=T, xlim = c(0,1), ylim = c(0,120))
# lines(x, u, col = "red")
# hist(values(cells.se)$cell_MF[values(cells.se)$cell_MF>0.6], breaks = 0.4*n_bins, prob=T, xlim = c(0,1), ylim = c(0,60))
# lines(x, m, col = "red")
# 




