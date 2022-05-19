library(SummarizedExperiment)
library(GenomicRanges)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
library(bsseq)
library(tidyverse)
setAutoBlockSize(size=1e9)
setwd(here::here())
n_bins <- 200


metadata <- fread("data/metadata/metadata_luo2017/NIHMS893063-supplement-Table_S2_csv.csv", skip = 1)
cell_counts <- metadata[, .(.N), by = .(`Neuron type`)] %>% arrange(desc(N)) 

dt_info <- cell_counts[`Neuron type` %in% c("hL2/3", "hL5-4", "hDL-1", "hPv-1", "hNdnf", "hSst-2")]
dt_info$CellClass <- c(rep("Excitatory", 3), rep("Inhibitory", 3))

dt_dens <- data.table()
for (i in 1:nrow(dt_info)) {  
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_human/subtype_", sub("/", "", dt_info$`Neuron type`[i]), 
                                                        "_", dt_info$N[i], "cells_qced"))
  dens <- calDensity(values(cells.se)$cell_MF, n_bins = n_bins)
  dt_dens <- rbind(dt_dens, data.table(Density = dens,
                                       CellClass = dt_info$CellClass[i],
                                       SubType = dt_info$`Neuron type`[i], 
                                       N_cell = dt_info$N[i],
                                       MedCovPct = metadata[`Neuron type`==dt_info$`Neuron type`[i], median(`Coverage (%)`)]))
  print(i)
}

dt_dens <- dt_dens %>%
  group_by(SubType) %>%
  mutate(CumDens = cumsum(Density.y))

write.table(dt_dens, "data/interim/estim_emiBetaPrior/eda/estimated_density_Luo2017_human_6subtypes.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")

dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Luo2017_human_6subtypes.txt") 
dt_dens %>% ggplot(aes(x = Density.x, y = CumDens)) + 
  geom_path(aes(linetype = CellClass, group = SubType, color = SubType)) +
  xlab("Value") +
  ylab("CDF of across-cell MF from subtypes in Luo2017 human data") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_Luo2017_human_6subtypes_coloredBySubtype.png", width = 7, height = 5)

dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Luo2017_human_6subtypes.txt") 
dt_dens %>% ggplot(aes(x = Density.x, y = CumDens)) + 
  geom_path(aes(linetype = CellClass, group = SubType, color = MedCovPct)) +
  viridis::scale_color_viridis(name = "Median coverage (%)", option = "D") +
  xlab("Value") +
  ylab("CDF of across-cell MF from subtypes in Luo2017 human data") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_Luo2017_human_6subtypes_coloredByCovPct.png", width = 7, height = 5)


dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Luo2017_human_6subtypes.txt") 
dt_dens %>% ggplot(aes(x = Density.x, y = CumDens)) + 
  geom_path(aes(linetype = CellClass, group = SubType, color = N_cell)) +
  viridis::scale_color_viridis(name = "N cells", option = "D") +
  xlab("Value") +
  ylab("CDF of across-cell MF from subtypes in Luo2017 human data") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_Luo2017_human_6subtypes_coloredByNcell.png", width = 7, height = 5)


dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Luo2017_human_6subtypes.txt") 
dt_dens %>% 
  pivot_wider(id_cols = Density.x, names_from = SubType, values_from = CumDens) %>%
  ggplot(aes(`hL2/3`, `hPv-1`)) +
  geom_point()
dt_dens %>%
  filter(SubType=="hL2/3") %>%
  ggplot(aes(Density.x, Density.y)) +
  geom_bar(stat = "identity")
dt_dens %>%
  filter(SubType=="hPv-1") %>%
  ggplot(aes(Density.x, Density.y)) +
  geom_bar(stat = "identity")


# ==== Effect of cell types on prior estimation (Investigate hL2/3 with N_cell = 873 and hPv-1 with N_cell = 175) ====
dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Luo2017_human_6subtypes.txt") 
x <- seq(from = 0, to = 1, by = 1/n_bins)
best_params_u_1 <- dt_dens %>%
  filter(SubType == "hL2/3") %>%
  select(1,2) %>%
  calBestParamFromDensity(population = "u", mc.cores = 16)
u1 <- 0.93 * pbeta(x, 1, 8) + 0.07 * pbeta(x, 1, 1720) 

best_params_m_1 <- dt_dens %>%
  filter(SubType == "hL2/3") %>%
  select(1,2) %>%
  calBestParamFromDensity(population = "m", mc.cores = 16)
m1 <- 0.97 * pbeta(x, 20, 1) + 0.03 * pbeta(x, 1650, 1)

best_params_u_2 <- dt_dens %>%
  filter(SubType == "hPv-1") %>%
  select(1,2) %>%
  calBestParamFromDensity(population = "u", mc.cores = 16)
u2 <- 0.83 * pbeta(x, 1, 6) + 0.17 * pbeta(x, 1, 1950) 

best_params_m_2 <- dt_dens %>%
  filter(SubType == "hPv-1") %>%
  select(1,2) %>%
  calBestParamFromDensity(population = "m", mc.cores = 16)
m2 <- 0.83 * pbeta(x, 7, 1) + 0.17 * pbeta(x, 1900, 1) 

data.frame(u1, u2) %>%
  ggplot(aes(u1, u2)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() +
  ggtitle("QQplot of fitted prior for unmethylated population from 2 subtypes") +
  xlab("Subtype hL2/3 with N_cell = 873") + ylab("Subtype hPv-1 with N_cell = 175")


data.frame(m1, m2) %>%
  ggplot(aes(m1, m2)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() +
  ggtitle("QQplot of fitted prior for methylated population from 2 subtypes") +
  xlab("Subtype hL2/3 with N_cell = 873") + ylab("Subtype hPv-1 with N_cell = 175")
