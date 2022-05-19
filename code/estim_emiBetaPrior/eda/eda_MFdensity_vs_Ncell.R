library(tidyverse)
library(SummarizedExperiment)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
setwd(here::here())

# ==== downsample N cells to half or quarter ====
cells.se <- loadHDF5SummarizedExperiment(dir = "data/processed/processed_luo2017_mice/subtype_mL23_649cells_qced")

cells_half.se <- cells.se[, 1:round(ncol(cells.se)/2)]
values(cells_half.se)$cell_cov <- rowSums(assays(cells_half.se)$Cov > 0, na.rm = T)
values(cells_half.se)$cell_meth <- rowSums(assays(cells_half.se)$M > 0, na.rm = T)
values(cells_half.se)$cell_MF <- values(cells_half.se)$cell_meth/values(cells_half.se)$cell_cov

cells_quart.se <- cells.se[, 1:round(ncol(cells.se)/4)]
values(cells_quart.se)$cell_cov <- rowSums(assays(cells_quart.se)$Cov > 0, na.rm = T)
values(cells_quart.se)$cell_meth <- rowSums(assays(cells_quart.se)$M > 0, na.rm = T)
values(cells_quart.se)$cell_MF <- values(cells_quart.se)$cell_meth/values(cells_quart.se)$cell_cov

png("plots/estim_emiBetaPrior/hist_MF_mL23_full_cells.png")
hist(values(cells.se)$cell_MF, breaks = 200, freq = F, main = "All cells", xlab = "MF")
dev.off()
png("plots/estim_emiBetaPrior/hist_MF_mL23_subsampled_half_cells.png")
hist(values(cells_half.se)$cell_MF, breaks = 200, freq = F, main = "Subsampled 1/2 cells", xlab = "MF")
dev.off()
png("plots/estim_emiBetaPrior/hist_MF_mL23_subsampled_quarter_cells.png")
hist(values(cells_quart.se)$cell_MF, breaks = 200, freq = F, main = "Subsampled 1/4 cells", xlab = "MF")
dev.off()

calBestParam(mf_vec = values(cells.se)$cell_MF, population = "u", mc.cores = 16)
# Fitted mixture distribution density for unmethylated population: 0.93 * dbeta(x, 1, 5) + 0.07 * dbeta(x, 1, 1860) 
# L1 difference between discrete empirical distribution and fitted mixture:0.24300550040935
calBestParam(mf_vec = values(cells_half.se)$cell_MF, population = "u", mc.cores = 16)
# Fitted mixture distribution density for unmethylated population: 0.9 * dbeta(x, 1, 4) + 0.1 * dbeta(x, 1, 1970) 
# L1 difference between discrete empirical distribution and fitted mixture:0.252664563157062
calBestParam(mf_vec = values(cells_quart.se)$cell_MF, population = "u", mc.cores = 16)
# Fitted mixture distribution density for unmethylated population: 0.85 * dbeta(x, 1, 3) + 0.15 * dbeta(x, 1, 1730) 
# L1 difference between discrete empirical distribution and fitted mixture:0.377466228214409

calBestParam(mf_vec = values(cells.se)$cell_MF, population = "u", param2_space = 1000, mc.cores = 16)
calBestParam(mf_vec = values(cells_half.se)$cell_MF, population = "u", param2_space = 1000, mc.cores = 16)
calBestParam(mf_vec = values(cells_quart.se)$cell_MF, population = "u", param2_space = 1000, mc.cores = 16)

calBestParam(mf_vec = values(cells.se)$cell_MF, population = "m", mc.cores = 16)
# Fitted mixture distribution density for methylated population: 0.97 * dbeta(x, 13, 1) + 0.03 * dbeta(x, 1670, 1) 
# L1 difference between discrete empirical distribution and fitted mixture:0.283478364323447
calBestParam(mf_vec = values(cells_half.se)$cell_MF, population = "m", mc.cores = 16)
# Fitted mixture distribution density for methylated population: 0.91 * dbeta(x, 8, 1) + 0.09 * dbeta(x, 1560, 1) 
# L1 difference between discrete empirical distribution and fitted mixture:0.366255849697933
calBestParam(mf_vec = values(cells_quart.se)$cell_MF, population = "m", mc.cores = 16)
# Fitted mixture distribution density for methylated population: 0.84 * dbeta(x, 5, 1) + 0.16 * dbeta(x, 1630, 1) 
# L1 difference between discrete empirical distribution and fitted mixture:0.443143237213531

dens_full <- calDensity(values(cells.se)$cell_MF, n_bins = 200)
dens_half <- calDensity(values(cells_half.se)$cell_MF, n_bins = 200)
ggplot() +
  geom_point(aes(c(0, cumsum(dens_full$y)), c(0, cumsum(dens_half$y)))) +
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  xlab("MF cumulative density calculated from all cells") +
  ylab("MF cumulative density calculated from subsampled half cells") +
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/qqplot_MF_mL23_fullVSsubsampled_cells.png", width = 6, height = 5)

