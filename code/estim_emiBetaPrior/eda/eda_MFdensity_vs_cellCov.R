library(tidyverse)
library(SummarizedExperiment)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
setwd(here::here())

summInfo <- function(df0, n_cell, dataset, cell_type){
  df <- df0 %>%
    as.data.frame() %>%
    group_by(cell_cov) %>%
    summarise(pct_MF_less0.5 = sum(cell_MF < 0.5)/length(cell_MF),
              pct_MF_greater0.5 = sum(cell_MF > 0.5)/length(cell_MF),
              mean_MF = mean(cell_MF), 
              sd_MF = sd(cell_MF),
              mean_MF_unmeth = mean(cell_MF[cell_MF < 0.5]),
              sd_MF_unmeth = sd(cell_MF[cell_MF < 0.5]),
              mean_MF_meth = mean(cell_MF[cell_MF > 0.5]),
              sd_MF_meth = sd(cell_MF[cell_MF > 0.5]),
              n_sites = length(cell_cov),
              n_sites_unmeth = length(cell_cov[cell_MF < 0.5]),
              n_sites_meth = length(cell_cov[cell_MF > 0.5]))
  df$n_cell <- n_cell
  df$dataset <- dataset
  df$cell_type <- cell_type
  return(df)
}

# ==== compare subtypes from different data source ====
## Exc subtypes in Liu2021
# cells_liu_exc_N1462 <- readRDS("data/processed/processed_liu2021/summary_subtype_DG_dg-all_1462cells_combined_qced.rds")
# cells_liu_exc_N1359 <- readRDS("data/processed/processed_liu2021/summary_subtype_CA1_Chrm3_1359cells_combined_qced.rds")
cells_liu_exc_N871.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample11F_190214_GSE135543_subtype_DG_dg-all_871cells_qced/"))
cells_liu_exc_N685.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample9H_190219_GSE135570_subtype_CA1_Chrm3_685cells_qced/"))
cells_liu_exc_N564.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample2E_180220_GSE131904_subtype_OLF-Exc_Bmpr1b_564cells_qced/"))
cells_liu_exc_N506.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample1A_180226_GSE130553_subtype_IT-L23_Cux1_506cells_qced/"))
## Inh subtypes in Liu2021
cells_liu_inh_N498.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample5H_181015_GSE134806_subtype_PAL-Inh_Meis2_498cells_qced/"))
cells_liu_inh_N439.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample1C_180208_GSE131836_subtype_OLF_Trpc4_439cells_qced/"))
cells_liu_inh_N369.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample5E_180925_GSE134653_subtype_MSN-D1_Hrh1_369cells_qced/"))
# cells_liu_inh_N353.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample4D_171219_GSE131406_subtype_MSN-D2_Slc24a2_353cells_qced/"))
## Exc subtypes in Luo2017
cells_luom_exc_N686.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_mice/subtype_mL6-2_686cells_qced/"))
cells_luom_exc_N649.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_mice/subtype_mL23_649cells_qced/"))
cells_luom_exc_N370.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_mice/subtype_mL4_370cells_qced/"))
## Inh subtypes in Luo2017
cells_luom_inh_N136.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_mice/subtype_mPv_136cells_qced/"))
cells_luom_inh_N123.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_mice/subtype_mSst-12_123cells_qced/"))

# cells_luoh_N873.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_human/subtype_hL23_873cells_qced/"))
# cells_luoh_N175.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_human/subtype_hPv-1_175cells_qced/"))
# cells_luoh_N162.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_human/subtype_hL5-4_162cells_qced/"))
# cells_luoh_N144.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_human/subtype_hDL-1_144cells_qced/"))

summ_info <- rbind(#summInfo(cells_liu_exc_N1462, 1462, "Liu2021", "Exc"),
                   #summInfo(cells_liu_exc_N1359, 1359, "Liu2021", "Exc"),
                   summInfo(values(cells_liu_exc_N871.se), 871, "Liu2021", "Exc"),
                   summInfo(values(cells_liu_exc_N685.se), 685, "Liu2021", "Exc"),
                   summInfo(values(cells_liu_exc_N564.se), 564, "Liu2021", "Exc"),
                   summInfo(values(cells_liu_exc_N506.se), 506, "Liu2021", "Exc"),
                   summInfo(values(cells_liu_inh_N498.se), 498, "Liu2021", "Inh"),
                   summInfo(values(cells_liu_inh_N439.se), 439, "Liu2021", "Inh"),
                   summInfo(values(cells_liu_inh_N369.se), 369, "Liu2021", "Inh"),
                   # summInfo(values(cells_liu_inh_N353.se), 353, "Liu2021", "Inh"),
                   summInfo(values(cells_luom_exc_N686.se), 686, "Luo2017m", "Exc"),
                   summInfo(values(cells_luom_exc_N649.se), 649, "Luo2017m", "Exc"),
                   summInfo(values(cells_luom_exc_N370.se), 370, "Luo2017m", "Exc"),
                   summInfo(values(cells_luom_inh_N136.se), 136, "Luo2017m", "Inh"),
                   summInfo(values(cells_luom_inh_N123.se), 123, "Luo2017m", "Inh")#,
                   # summInfo(cells_luoh_N873.se, "Luo2017h"),
                   # summInfo(cells_luoh_N175.se, "Luo2017h"),
                   # summInfo(cells_luoh_N162.se, "Luo2017h"),
                   # summInfo(cells_luoh_N144.se, "Luo2017h")
                   )
summ_info$dataset <- factor(summ_info$dataset, levels = c("Liu2021", "Luo2017m", "Luo2017h"))

## mean MF vs cell cov
summ_info %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(x = cell_cov/n_cell, y = mean_MF, color = n_cell, shape = dataset)) + 
  geom_errorbar(aes(ymin = sapply(mean_MF - sd_MF, function(x) max(0,x)), 
                    ymax = sapply(mean_MF + sd_MF, function(x) min(1,x))),
                alpha = 0.1, size = 1, color = "grey") +
  geom_smooth(aes(group = n_cell, linetype = dataset), method = "loess", span = 0.5, size = 0.7, se = F) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  theme_bw() +
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_meanMF_vs_cellCov_coloredByNcells.png", width = 7, height = 5)
summ_info %>%
  filter(n_sites >= 10000) %>%
  group_by(n_cell, dataset) %>%
  mutate(cdf_n_sites = cumsum(n_sites)/sum(n_sites)) %>%
  ggplot(aes(cell_cov/n_cell, cdf_n_sites, color = n_cell, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), method = "loess", span = 0.5, size = 0.7, se = F) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  theme_bw() +
  xlab("Across-cell coverage / # cells") + ylab("CDF of # sites")
ggsave("plots/estim_emiBetaPrior/pointSmooth_NsitesCDF_vs_cellCov_coloredByNcells.png", width = 7, height = 5)

## percent of MF<0.5
summ_info %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, pct_MF_less0.5, color = n_cell, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), method = "loess", span = 0.5, size = 0.7, se = F) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  theme_bw() +
  ggtitle("Low methylated sites (MF < 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Percentage of sites with MF<0.5")
ggsave("plots/estim_emiBetaPrior/pointSmooth_pctMFless0.5_vs_cellCov_coloredByNcells.png", width = 7, height = 5)

## percent of MF>0.5
summ_info %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, pct_MF_greater0.5, color = n_cell, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), method = "loess", span = 0.5, size = 0.7, se = F) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  theme_bw() + 
  ggtitle("High methylated sites (MF > 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Percentage of sites with MF>0.5")
ggsave("plots/estim_emiBetaPrior/pointSmooth_pctMFgreater0.5_vs_cellCov_coloredByNcells.png", width = 7, height = 5)


## mean MF vs cell cov in unmeth subpop colored by number of cells
summ_info %>%
  filter(n_sites_unmeth >= 10000) %>%
  # filter(dataset == "Liu2021") %>%
  ggplot(aes(x = cell_cov/n_cell, y = mean_MF_unmeth, group = n_cell, shape = dataset)) + 
  geom_errorbar(aes(ymin = sapply(mean_MF_unmeth - sd_MF_unmeth, function(x) max(0,x)), 
                    ymax = sapply(mean_MF_unmeth + sd_MF_unmeth, function(x) min(1,x))),
                alpha = 0.2, size = 1, color = "grey") +
  geom_smooth(aes(color = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  # geom_line(aes(color = cell_type, linetype = dataset)) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  theme_bw() + ylim(0,0.5) +
  ggtitle("Low methylated sites (MF < 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_meanMF_vs_cellCov_unmethPop_coloredByNcells.png", width = 7, height = 5)
summ_info %>%
  filter(n_sites_unmeth >= 10000) %>%
  group_by(n_cell, dataset) %>%
  mutate(cdf_n_sites_unmeth = cumsum(n_sites_unmeth)/sum(n_sites_unmeth)) %>%
  ggplot(aes(cell_cov/n_cell, cdf_n_sites_unmeth, color = n_cell, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  theme_bw() + 
  ggtitle("Low methylated sites (MF < 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("CDF of # sites")
ggsave("plots/estim_emiBetaPrior/pointSmooth_NsitesCDF_vs_cellCov_unmethPop_coloredByNcells.png", width = 7, height = 5)
summ_info %>%
  filter(n_sites_unmeth >= 10000) %>%
  group_by(n_cell, dataset) %>%
  ggplot(aes(cell_cov/n_cell, n_sites_unmeth, color = n_cell, shape = dataset)) + 
  geom_path(aes(group = n_cell, linetype = dataset), size = 0.7) +
  geom_abline(intercept = 10000, slope = 0, color = "grey", linetype = "dashed") + 
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  theme_bw() + 
  ggtitle("Low methylated sites (MF < 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("# sites")
ggsave("plots/estim_emiBetaPrior/line_Nsites_vs_cellCov_unmethPop_coloredByNcells.png", width = 7, height = 5)

## mean MF vs cell cov in meth subpop colored by number of cells
summ_info %>%
  filter(n_sites_meth >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, mean_MF_meth, color = n_cell)) + 
  geom_errorbar(aes(ymin = sapply(mean_MF_meth - sd_MF_meth, function(x) max(0,x)), 
                    ymax = sapply(mean_MF_meth + sd_MF_meth, function(x) min(1,x))),
                alpha = 0.1, size = 1, color = "grey") +
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  theme_bw() + ylim(0.5,1) +
  ggtitle("High methylated sites (MF > 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_meanMF_vs_cellCov_methPop_coloredByNcells.png", width = 7, height = 5)
summ_info %>%
  filter(n_sites_meth >= 10000) %>%
  group_by(n_cell, dataset) %>%
  mutate(cdf_n_sites_meth = cumsum(n_sites_meth)/sum(n_sites_meth)) %>%
  ggplot(aes(cell_cov/n_cell, cdf_n_sites_meth, color = n_cell)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  theme_bw() +
  ggtitle("High methylated sites (MF > 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("CDF of # sites")
ggsave("plots/estim_emiBetaPrior/pointSmooth_NsitesCDF_vs_cellCov_methPop_coloredByNcells.png", width = 7, height = 5)
summ_info %>%
  filter(n_sites_meth >= 10000) %>%
  group_by(n_cell, dataset) %>%
  ggplot(aes(cell_cov/n_cell, n_sites_meth, color = n_cell, shape = dataset)) + 
  geom_path(aes(group = n_cell, linetype = dataset), size = 0.7) +
  geom_abline(intercept = 10000, slope = 0, color = "grey", linetype = "dashed") + 
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  theme_bw() + 
  ggtitle("High methylated sites (MF > 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("# sites")
ggsave("plots/estim_emiBetaPrior/line_Nsites_vs_cellCov_methPop_coloredByNcells.png", width = 7, height = 5)

## mean MF vs cell cov in unmeth subpop colored by cell type
summ_info %>%
  filter(n_sites >= 10000) %>%
  # filter(dataset == "Liu2021") %>%
  ggplot(aes(x = cell_cov/n_cell, y = mean_MF_unmeth, group = n_cell)) + 
  geom_errorbar(aes(ymin = sapply(mean_MF_unmeth - sd_MF_unmeth, function(x) max(0,x)), 
                    ymax = sapply(mean_MF_unmeth + sd_MF_unmeth, function(x) min(1,x))),
                alpha = 0.2, size = 1, color = "grey") +
  geom_smooth(aes(color = cell_type, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  theme_bw() + ylim(0,0.5) +
  ggtitle("Low methylated sites (MF < 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_meanMF_vs_cellCov_unmethPop_coloredByCellType.png", width = 7, height = 5)


## mean MF vs cell cov in meth subpop colored by cell type
summ_info %>%
  filter(n_sites >= 10000) %>%
  # filter(dataset == "Liu2021") %>%
  ggplot(aes(x = cell_cov/n_cell, y = mean_MF_meth, group = n_cell)) + 
  geom_errorbar(aes(ymin = sapply(mean_MF_meth - sd_MF_meth, function(x) max(0,x)), 
                    ymax = sapply(mean_MF_meth + sd_MF_meth, function(x) min(1,x))),
                alpha = 0.2, size = 1, color = "grey") +
  geom_smooth(aes(color = cell_type, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  theme_bw() + ylim(0.5, 1) +
  ggtitle("High methylated sites (MF > 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_meanMF_vs_cellCov_methPop_coloredByCellType.png", width = 7, height = 5)


# ==== compare subsamples of subtype DG dg-all in Liu2021 data ====
cells_liu_N871.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample11F_190214_GSE135543_subtype_DG_dg-all_871cells_qced/"))
half_liu_N871.se <- cells_liu_N871.se[,1:floor(ncol(cells_liu_N871.se)/2)]
quart_liu_N871.se <- cells_liu_N871.se[,1:floor(ncol(cells_liu_N871.se)/4)]

values(half_liu_N871.se)$cell_cov <- rowSums(assays(half_liu_N871.se)$Cov > 0, na.rm = T)
values(half_liu_N871.se)$cell_meth <- rowSums(assays(half_liu_N871.se)$M > 0, na.rm = T)
values(half_liu_N871.se)$cell_MF <- values(half_liu_N871.se)$cell_meth/values(half_liu_N871.se)$cell_cov
half_liu_N871.se <- half_liu_N871.se[values(half_liu_N871.se)$cell_cov>=5,]

values(quart_liu_N871.se)$cell_cov <- rowSums(assays(quart_liu_N871.se)$Cov > 0, na.rm = T)
values(quart_liu_N871.se)$cell_meth <- rowSums(assays(quart_liu_N871.se)$M > 0, na.rm = T)
values(quart_liu_N871.se)$cell_MF <- values(quart_liu_N871.se)$cell_meth/values(quart_liu_N871.se)$cell_cov
quart_liu_N871.se <- quart_liu_N871.se[values(quart_liu_N871.se)$cell_cov>=5,]

summ_info2 <- rbind(summInfo(values(cells_liu_N871.se), "Full"),
                    summInfo(values(half_liu_N871.se), "Half"),
                    summInfo(values(quart_liu_N871.se), "Quarter")
                   )

## mean MF vs cell cov
summ_info2 %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, mean_MF, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() +
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Liu2021_dg-all_subsamples_meanMF_vs_cellCov_coloredByNcells.png", width = 7, height = 5)

## percent of MF<0.5
summ_info2 %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, pct_MF_less0.5, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() +
  xlab("Across-cell coverage / # cells") + ylab("Percentage of sites with MF<0.5")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Liu2021_dg-all_subsamples_pctMFless0.5_vs_cellCov_coloredByNcells.png", width = 7, height = 5)

## percent of MF>0.5
summ_info2 %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, pct_MF_greater0.5, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() + 
  xlab("Across-cell coverage / # cells") + ylab("Percentage of sites with MF>0.5")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Liu2021_dg-all_subsamples_pctMFgreater0.5_vs_cellCov_coloredByNcells.png", width = 7, height = 5)


## mean MF vs cell cov in unmeth subpop
summ_info2 %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, mean_MF_unmeth, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() + ylim(0, 0.5) +
  ggtitle("Low methylated sites (MF < 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Liu2021_dg-all_subsamples_meanMF_vs_cellCov_unmethPop_coloredByNcells.png", width = 7, height = 5)
summ_info2 %>%
  filter(n_sites >= 10000) %>%
  group_by(n_cell, dataset) %>%
  mutate(cdf_n_sites_unmeth = cumsum(n_sites_unmeth)/sum(n_sites_unmeth)) %>%
  ggplot(aes(cell_cov/n_cell, cdf_n_sites_unmeth, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(3,2,1), name = "Dataset source") +
  theme_bw() +
  ggtitle("Low methylated sites (MF < 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("CDF of # sites")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Liu2021_dg-all_subsamples_NsitesCDF_vs_cellCov_unmethPop_coloredByNcells.png", width = 7, height = 5)

## mean MF vs cell cov in meth subpop
summ_info2 %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, mean_MF_meth, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() + ylim(0.5, 1) +
  ggtitle("High methylated sites (MF > 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Liu2021_dg-all_subsamples_meanMF_vs_cellCov_methPop_coloredByNcells.png", width = 7, height = 5)
summ_info2 %>%
  filter(n_sites >= 10000) %>%
  group_by(n_cell, dataset) %>%
  mutate(cdf_n_sites_meth = cumsum(n_sites_meth)/sum(n_sites_meth)) %>%
  ggplot(aes(cell_cov/n_cell, cdf_n_sites_meth, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() +
  ggtitle("High methylated sites (MF > 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("CDF of # sites")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Liu2021_dg-all_subsamples_NsitesCDF_vs_cellCov_methPop_coloredByNcells.png", width = 7, height = 5)


# ==== compare subsamples of subtype mL2/3 in Luo2017 mice data ====
cells_luom_N649.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_mice/subtype_mL23_649cells_qced/"))
half_luom_N649.se <- cells_luom_N649.se[,1:floor(ncol(cells_luom_N649.se)/2)]
quart_luom_N649.se <- cells_luom_N649.se[,1:floor(ncol(cells_luom_N649.se)/4)]

values(half_luom_N649.se)$cell_cov <- rowSums(assays(half_luom_N649.se)$Cov > 0, na.rm = T)
values(half_luom_N649.se)$cell_meth <- rowSums(assays(half_luom_N649.se)$M > 0, na.rm = T)
values(half_luom_N649.se)$cell_MF <- values(half_luom_N649.se)$cell_meth/values(half_luom_N649.se)$cell_cov
half_luom_N649.se <- half_luom_N649.se[values(half_luom_N649.se)$cell_cov>=5,]

values(quart_luom_N649.se)$cell_cov <- rowSums(assays(quart_luom_N649.se)$Cov > 0, na.rm = T)
values(quart_luom_N649.se)$cell_meth <- rowSums(assays(quart_luom_N649.se)$M > 0, na.rm = T)
values(quart_luom_N649.se)$cell_MF <- values(quart_luom_N649.se)$cell_meth/values(quart_luom_N649.se)$cell_cov
quart_luom_N649.se <- quart_luom_N649.se[values(quart_luom_N649.se)$cell_cov>=5,]

summ_info3 <- rbind(summInfo(values(cells_luom_N649.se), "Full"),
                    summInfo(values(half_luom_N649.se), "Half"),
                    summInfo(values(quart_luom_N649.se), "Quarter")
)

## mean MF vs cell cov
summ_info3 %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, mean_MF, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() +
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Luo2017_mL23_subsamples_meanMF_vs_cellCov_coloredByNcells.png", width = 7, height = 5)

## percent of MF<0.5
summ_info3 %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, pct_MF_less0.5, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() +
  xlab("Across-cell coverage / # cells") + ylab("Percentage of sites with MF<0.5")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Luo2017_mL23_subsamples_pctMFless0.5_vs_cellCov_coloredByNcells.png", width = 7, height = 5)

## percent of MF>0.5
summ_info3 %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, pct_MF_greater0.5, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() +
  xlab("Across-cell coverage / # cells") + ylab("Percentage of sites with MF>0.5")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Luo2017_mL23_subsamples_pctMFgreater0.5_vs_cellCov_coloredByNcells.png", width = 7, height = 5)


## mean MF vs cell cov in unmeth subpop
summ_info3 %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, mean_MF_unmeth, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() + ylim(0, 0.5) +
  ggtitle("Low methylated sites (MF < 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Luo2017_mL23_subsamples_meanMF_vs_cellCov_unmethPop_coloredByNcells.png", width = 7, height = 5)
summ_info3 %>%
  filter(n_sites >= 10000) %>%
  group_by(n_cell, dataset) %>%
  mutate(cdf_n_sites_unmeth = cumsum(n_sites_unmeth)/sum(n_sites_unmeth)) %>%
  ggplot(aes(cell_cov/n_cell, cdf_n_sites_unmeth, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() +
  ggtitle("Low methylated sites (MF < 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("CDF of # sites")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Luo2017_mL23_subsamples_NsitesCDF_vs_cellCov_unmethPop_coloredByNcells.png", width = 7, height = 5)

## mean MF vs cell cov in meth subpop
summ_info3 %>%
  filter(n_sites >= 10000) %>%
  ggplot(aes(cell_cov/n_cell, mean_MF_meth, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() +  ylim(0.5, 1) +
  ggtitle("High methylated sites (MF > 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Luo2017_mL23_subsamples_meanMF_vs_cellCov_methPop_coloredByNcells.png", width = 7, height = 5)
summ_info3 %>%
  filter(n_sites >= 10000) %>%
  group_by(n_cell, dataset) %>%
  mutate(cdf_n_sites_meth = cumsum(n_sites_meth)/sum(n_sites_meth)) %>%
  ggplot(aes(cell_cov/n_cell, cdf_n_sites_meth, shape = dataset)) + 
  geom_smooth(aes(group = n_cell, linetype = dataset), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  scale_linetype_manual(values=c(1,2,3), name = "Dataset source") +
  theme_bw() +
  ggtitle("High methylated sites (MF > 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("CDF of # sites")
ggsave("plots/estim_emiBetaPrior/pointSmooth_Luo2017_mL23_subsamples_NsitesCDF_vs_cellCov_methPop_coloredByNcells.png", width = 7, height = 5)
