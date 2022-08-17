library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")

# ==== combining CA1 Chrm3 cells from sample 9H_190219 and 9H_190212 ====
se_obj1 <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample9H_190219_GSE135570_subtype_CA1_Chrm3_685cells/"))
se_obj2 <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample9H_190212_GSE135447_subtype_CA1_Chrm3_674cells/"))

# compute across-cell coverage and methylated counts
values(se_obj1)$cell_cov <- rowSums(assays(se_obj1)$Cov > 0, na.rm = T)
values(se_obj1)$cell_meth <- rowSums(assays(se_obj1)$M > 0, na.rm = T)
values(se_obj2)$cell_cov <- rowSums(assays(se_obj2)$Cov > 0, na.rm = T)
values(se_obj2)$cell_meth <- rowSums(assays(se_obj2)$M > 0, na.rm = T)

df_obj1 <- as.data.table(granges(se_obj1))
df_obj2 <- as.data.table(granges(se_obj2))
df_objs <- merge.data.table(df_obj1, df_obj2, by = c("seqnames","start","end","width","strand"), 
                            suffixes = c("1", "2"), all = T) %>%
  mutate(cell_cov1 = replace_na(cell_cov1, 0), 
         cell_meth1 = replace_na(cell_meth1, 0),
         cell_cov2 = replace_na(cell_cov2, 0),
         cell_meth2 = replace_na(cell_meth2, 0)) %>% ## fill 0's to missing elements
  mutate(cell_cov = cell_cov1 + cell_cov2, cell_meth = cell_meth1 + cell_meth2)%>%
  filter(cell_cov >= 5) %>%
  mutate(cell_MF = cell_meth / cell_cov) 

saveRDS(df_objs %>%
          select(-c("cell_cov1", "cell_meth1", "cell_cov2", "cell_meth2")), 
        file = paste0("data/processed/processed_liu2021/summary_subtype_CA1_Chrm3_", ncol(se_obj1)+ncol(se_obj2), "cells_combined_qced.rds"))

df_objs %>%
  mutate(cell_MF1_bins = cut(cell_meth1 / cell_cov1, breaks = 50)) %>%
  ggplot(aes(x = cell_MF1_bins, y = cell_meth2/cell_cov2)) +
  geom_boxplot()



# ==== combining DG dg-all cells from sample 11F_190214 and 8J_190716 ====
se_obj1 <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample11F_190214_GSE135543_subtype_DG_dg-all_871cells/"))
se_obj2 <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample8J_190716_GSE158154_subtype_DG_dg-all_591cells/"))

# compute across-cell coverage and methylated counts
values(se_obj1)$cell_cov <- rowSums(assays(se_obj1)$Cov > 0, na.rm = T)
values(se_obj1)$cell_meth <- rowSums(assays(se_obj1)$M > 0, na.rm = T)
values(se_obj2)$cell_cov <- rowSums(assays(se_obj2)$Cov > 0, na.rm = T)
values(se_obj2)$cell_meth <- rowSums(assays(se_obj2)$M > 0, na.rm = T)

df_obj1 <- as.data.table(granges(se_obj1))
df_obj2 <- as.data.table(granges(se_obj2))
df_objs <- merge.data.table(df_obj1, df_obj2, by = c("seqnames","start","end","width","strand"), 
                            suffixes = c("1", "2"), all = T) %>%
  mutate(cell_cov1 = replace_na(cell_cov1, 0), 
         cell_meth1 = replace_na(cell_meth1, 0),
         cell_cov2 = replace_na(cell_cov2, 0),
         cell_meth2 = replace_na(cell_meth2, 0)) %>% ## fill 0's to missing elements
  mutate(cell_cov = cell_cov1 + cell_cov2, cell_meth = cell_meth1 + cell_meth2) %>%
  select(-c("cell_cov1", "cell_meth1", "cell_cov2", "cell_meth2")) %>%
  filter(cell_cov >= 5) %>%
  mutate(cell_MF = cell_meth / cell_cov)

saveRDS(df_objs, file = paste0("data/processed/processed_liu2021/summary_subtype_DG_dg-all_", ncol(se_obj1)+ncol(se_obj2), "cells_combined_qced.rds"))


# ==== before & after combining ====

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

cells_liu_exc_N1462 <- readRDS("data/processed/processed_liu2021/summary_subtype_DG_dg-all_1462cells_combined_qced.rds")
cells_liu_exc_N871.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample11F_190214_GSE135543_subtype_DG_dg-all_871cells_qced/"))
cells_liu_exc_N591.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample8J_190716_GSE158154_subtype_DG_dg-all_591cells_qced//"))

cells_liu_exc_N1359 <- readRDS("data/processed/processed_liu2021/summary_subtype_CA1_Chrm3_1359cells_combined_qced.rds")
cells_liu_exc_N685.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample9H_190219_GSE135570_subtype_CA1_Chrm3_685cells_qced/"))
cells_liu_exc_N674.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample9H_190212_GSE135447_subtype_CA1_Chrm3_674cells_qced//"))


summ_info0 <- rbind(summInfo(cells_liu_exc_N1462, 1462, "Liu2021", "DG dg-all"),
                   summInfo(cells_liu_exc_N1359, 1359, "Liu2021", "CA1 Chrm3"),
                   summInfo(values(cells_liu_exc_N871.se), 871, "Liu2021", "DG dg-all"),
                   summInfo(values(cells_liu_exc_N591.se), 591, "Liu2021", "DG dg-all"),
                   summInfo(values(cells_liu_exc_N685.se), 685, "Liu2021", "CA1 Chrm3"),
                   summInfo(values(cells_liu_exc_N674.se), 674, "Liu2021", "CA1 Chrm3")
                   )

summ_info0 %>%
  filter(n_sites_unmeth >= 10000, cell_cov >= 5) %>%
  ggplot(aes(x = cell_cov/n_cell, y = mean_MF_unmeth, group = n_cell)) + 
  geom_errorbar(aes(ymin = sapply(mean_MF_unmeth - sd_MF_unmeth, function(x) max(0,x)), 
                    ymax = sapply(mean_MF_unmeth + sd_MF_unmeth, function(x) min(1,x))),
                alpha = 0.2, size = 1, color = "grey") +
  geom_smooth(aes(color = n_cell, linetype = cell_type), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  # geom_line(aes(color = cell_type, linetype = dataset)) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  theme_bw() + ylim(0,0.5) + xlim(0, 0.12) +
  ggtitle("Low methylated sites (MF < 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_sampleCombining_beforeVSafter_meanMF_vs_cellCov_unmethPop.png", width = 7, height = 5)


summ_info0 %>%
  filter(n_sites_meth >= 10000, cell_cov >= 5) %>%
  ggplot(aes(x = cell_cov/n_cell, y = mean_MF_meth, group = n_cell)) + 
  geom_errorbar(aes(ymin = sapply(mean_MF_meth - sd_MF_meth, function(x) max(0,x)), 
                    ymax = sapply(mean_MF_meth + sd_MF_meth, function(x) min(1,x))),
                alpha = 0.2, size = 1, color = "grey") +
  geom_smooth(aes(color = n_cell, linetype = cell_type), size = 0.7, method = "loess", span = 0.5, se = FALSE) +
  # geom_line(aes(color = cell_type, linetype = dataset)) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  scale_linetype_manual(values=c(1,2), name = "Dataset source") +
  theme_bw() + ylim(0.5, 1) +
  ggtitle("High methylated sites (MF > 0.5)") + 
  xlab("Across-cell coverage / # cells") + ylab("Average MF given coverage")
ggsave("plots/estim_emiBetaPrior/pointSmooth_sampleCombining_beforeVSafter_meanMF_vs_cellCov_methPop.png", width = 7, height = 5)
