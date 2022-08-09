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
metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv")
cell_counts <- metadata[, .(.N), by = .(Sample, CellClass, SubType, GEO_accession)]
metadata[, .(.N), by = .(CellClass, SubType)] %>%  
  arrange(desc(N)) 

dt_info_Exc <- cell_counts %>%
  filter(CellClass=="Exc") %>%  
  arrange(desc(N)) %>%
  filter(!is.na(SubType) & N >= 500 & N <= 871)
dt_info_Exc <- dt_info_Exc[-(3:4),] # remove repeating subtypes
dt_info_Exc
#        Sample CellClass        SubType GEO_accession   N
# 1: 11F_190214       Exc      DG dg-all     GSE135543 871
# 2:  9H_190219       Exc      CA1 Chrm3     GSE135570 685
# 3:  2E_180220       Exc OLF-Exc Bmpr1b     GSE131904 564
# 4:  1A_180226       Exc    IT-L23 Cux1     GSE130553 506

dt_info_Inh <- cell_counts %>%
  filter(CellClass=="Inh") %>%  
  arrange(desc(N)) %>%
  filter(!is.na(SubType) & N >= 350)
dt_info_Inh <- dt_info_Inh[-c(3,5),] # remove repeating subtypes
dt_info_Inh
#       Sample CellClass        SubType GEO_accession   N
# 1: 5H_181015       Inh  PAL-Inh Meis2     GSE134806 498
# 2: 1C_180208       Inh      OLF Trpc4     GSE131836 439
# 3: 5E_180925       Inh    MSN-D1 Hrh1     GSE134653 369
# 4: 5E_180925       Inh MSN-D2 Slc24a2     GSE134653 352

dt_info <- rbind(dt_info_Exc, dt_info_Inh)


# ==== Estimate MF density ====
dt_dens <- data.table()
for (i in 1:nrow(dt_info)) {
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample", dt_info$Sample[i],
                                                        "_", dt_info$GEO_accession[i], "_subtype_", sub(" ", "_", dt_info$SubType[i]), 
                                                        "_", dt_info$N[i], "cells_qced"))
  dens <- calDensity(values(cells.se)$cell_MF, n_bins = n_bins)
  params_u <- calBestParamFromDensity(dens, population = "u", param2_space = seq(1000,2000,10), mc.cores = 16)
  params_m <- calBestParamFromDensity(dens, population = "m", param2_space = seq(1000,2000,10), mc.cores = 16)
  
  dt_dens <- rbind(dt_dens, data.table(Density = dens, 
                                       Sample = dt_info$Sample[i], 
                                       GEO_accession = dt_info$GEO_accession[i], 
                                       CellClass = dt_info$CellClass[i],
                                       SubType = dt_info$SubType[i], 
                                       N_cell = dt_info$N[i],
                                       MedFinalReads = metadata[Sample==dt_info$Sample[i], median(FinalReads)/1e6],
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
  group_by(Sample, SubType) %>%
  mutate(CumDens = cumsum(Density.y))

write.table(dt_dens, "data/interim/estim_emiBetaPrior/eda/estimated_density_Liu2021_8subtypes.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")

# ==== Plot cumulative density ====
dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Liu2021_8subtypes.txt") 
dt_dens %>% ggplot(aes(x = Density.x, y = CumDens)) + 
  geom_path(aes(linetype = CellClass, group = SubType, color = MedFinalReads)) +
  # scale_color_gradient(high = "red", low = "blue") +
  viridis::scale_color_viridis(name = "Median cell\ntotal reads\n/10^6", option = "D") +
  xlab("Value") +
  ylab("CDF of fractional methylation from subtypes in Liu2021 data") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_Liu2021_8subtypes_coloredByMedReads.png", width = 7, height = 5)

dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Liu2021_8subtypes.txt") 
dt_dens %>% ggplot(aes(x = Density.x, y = CumDens)) + 
  geom_path(aes(linetype = CellClass, color = SubType)) +
  # scale_color_gradient(high = "red", low = "blue") +
  xlab("Value") +
  ylab("CDF of fractional methylation from subtypes in Liu2021 data") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_Liu2021_8subtypes_coloredBySubtype.png", width = 7, height = 5)

dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Liu2021_8subtypes.txt") 
dt_dens %>% ggplot(aes(x = Density.x, y = CumDens)) + 
  geom_path(aes(linetype = CellClass, group = SubType, color = N_cell)) +
  # scale_color_gradient(high = "red", low = "blue") +
  viridis::scale_color_viridis(name = "N cells", option = "D") +
  xlab("Value") +
  ylab("CDF of fractional methylation from subtypes in Liu2021 data") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_Liu2021_8subtypes_coloredByNcell.png", width = 7, height = 5)

dt_dens <- fread("data/interim/estim_emiBetaPrior/eda/estimated_density_Liu2021_8subtypes.txt") 
dt_dens %>%
  filter(SubType=="MSN-D2 Slc24a2") %>%
  ggplot(aes(Density.x, Density.y)) +
  geom_bar(stat = "identity")
dt_dens %>%
  filter(SubType=="DG dg-all") %>%
  ggplot(aes(Density.x, Density.y)) +
  geom_bar(stat = "identity")
