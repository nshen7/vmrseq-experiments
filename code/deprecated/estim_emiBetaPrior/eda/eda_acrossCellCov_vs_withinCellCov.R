library(SummarizedExperiment)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
setwd(here::here())

# ==== helper functions ====
my_fun <- function(block)
  apply(block, 1, FUN = function(r){mean(r[r>0])}) # calculate the mean within-cell coverage for each site

my_rowMeans <- function(x) {
  grid <- rowAutoGrid(x = x, nrow = 1e6)
  block_level_stat <- blockApply(x, FUN = my_fun, grid = grid)
  unlist(block_level_stat)
}

# ==== across-cell coverage vs. within-cell coverage (liu2021 DG_dg-all 871cells) ====
cells_liu_N871.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample11F_190214_GSE135543_subtype_DG_dg-all_871cells_qced/"))

values(cells_liu_N871.se)$avg_cov_within_cell <- 
  my_rowMeans(assays(cells_liu_N871.se)$Cov)


values(cells_liu_N871.se) %>%
  as.data.frame() %>%
  group_by(cell_cov) %>%
  summarise(mean = mean(avg_cov_within_cell),
            sd = sd(avg_cov_within_cell)) %>%
  ggplot(aes(x = cell_cov, y = mean)) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width=.2,
                position = position_dodge(.9), color = "grey") +
  geom_point(size = 0.5) +
  theme_bw() +
  ylab("Average within-cell coverage") +
  xlab("Across-cell coverage") +
  ggtitle("Liu2021 DG_dg-all 871cells")
ggsave("plots/estim_emiBetaPrior/eda/point_Liu2021_dg-all_within_vs_across_cellCov.png")

# ==== across-cell coverage vs. within-cell coverage (luo2017 mL4 370cells) ====
cells_luom_N370.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_mice/subtype_mL4_370cells_qced/"))

values(cells_luom_N370.se)$avg_cov_within_cell <- 
  my_rowMeans(assays(cells_luom_N370.se)$Cov)

values(cells_luom_N370.se) %>%
  as.data.frame() %>%
  group_by(cell_cov) %>%
  summarise(mean = mean(avg_cov_within_cell),
            sd = sd(avg_cov_within_cell)) %>%
  ggplot(aes(x = cell_cov, y = mean)) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width=.2,
                position = position_dodge(.9), color = "grey") +
  geom_point(size = 0.5) +
  theme_bw() +
  ylab("Average within-cell coverage") +
  xlab("Across-cell coverage") +
  ggtitle("luo2017 mL4 370 cells")
ggsave("plots/estim_emiBetaPrior/eda/point_luo2017_mL4_within_vs_across_cellCov.png")
