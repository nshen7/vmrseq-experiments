library(data.table)
library(readxl)
library(tidyverse)
library(SummarizedExperiment)
library(GenomicRanges)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
library(bsseq)
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")

metadata <- fread("../../DXM_extend_chr1/data/metadata/sample_info_processed.csv")
table(metadata$Neuron_type3)
metadata[, median(`Coverage(%)`), by = Neuron_type3]

metadata %>%
  ggplot(aes(`Coverage(%)`, `mCG/CG`, color = Neuron_type1)) +
  geom_point(size = 0.5, alpha = 0.5) + 
  geom_smooth(method = "lm") +
  theme_bw()

metadata %>%
  ggplot(aes(Neuron_type3, `Coverage(%)`)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Subtypes") 
ggsave("plots/eda_luo2017_mice/boxplot_readCoverage_vs_subtypes.png", width = 7, height = 5)
### all cell types here have comparable read coverage

# ==== total coverage per cell vs. MF ====
total_cov_mL23 <- colSums(assays(cells_mL23.se)$Cov > 0, na.rm = T)
total_meth_mL23 <- colSums(assays(cells_mL23.se)$M/assays(cells_mL23.se)$Cov >= 0.5, na.rm = T) 
total_cov_mL62 <- colSums(assays(cells_mL62.se)$Cov > 0, na.rm = T)
total_meth_mL62 <- colSums(assays(cells_mL62.se)$M/assays(cells_mL62.se)$Cov >= 0.5, na.rm = T) 

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
data.frame(total_cov = c(total_cov_mL23, total_cov_mL62),
           total_meth = c(total_meth_mL23, total_meth_mL62),
           cov_level = c(rep("mL2/3", length(total_cov_mL23)), 
                         rep("mL6-2", length(total_cov_mL62))
           )) %>%
  mutate(total_MF = total_meth / total_cov) %>%
  ggplot(aes(total_cov, total_MF, color = cov_level)) +
  geom_point(size = 1) + 
  geom_smooth(aes(group = cov_level), method = "lm") + 
  scale_color_manual(values = cbp1, name="Sample") +
  theme_bw() +
  xlab("Total covered sites in cell") +
  ylab("mCG/CG") + ylim(0.6,0.9)
ggsave("plots/eda_luo2017_mice/point_MFofSites_vs_totalCoveredSites_2Subtypes.png", width = 7, height = 5)  

# ==== pre-QC mL2/3 649 cells ====
cells_mL23.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_mice/subtype_mL23_649cells/"))
values(cells_mL23.se)$cell_cov <- rowSums(assays(cells_mL23.se)$Cov > 0, na.rm = T)
values(cells_mL23.se)$cell_meth <- rowSums(assays(cells_mL23.se)$M/assays(cells_mL23.se)$Cov >= 0.5, na.rm = T) 
values(cells_mL23.se)$cell_MF <- values(cells_mL23.se)$cell_meth/values(cells_mL23.se)$cell_cov

# hist(values(cells.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_mL23.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_mL23.se)$cell_cov, xlim = c(0,100), breaks = 300, freq = F)
quantile(values(cells_mL23.se)$cell_cov)
quantile(values(cells_mL23.se)$cell_cov, probs = 0.05)
values(cells_mL23.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.05)))
  ) %>%
  ggplot(aes(cell_cov_quantile, cell_MF)) +
  geom_boxplot() +
  theme_bw() + 
  ggtitle("Subtype mL2/3 - 649 cells") +
  xlab("Cell coverage quantile") + ylab("Cell MF") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_luo2017_mice/boxplot_cellMF_vs_QcellCov_subtype_mL23_649cells.png", width = 8, height = 5)

values(cells_mL23.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  slice_sample(prop = 0.1) %>% # suubset a random 10% of the sites to quickly plot
  filter(cell_cov_percent <= 0.15) %>%
  ggplot(aes(cell_cov_percent, cell_MF)) +
  stat_density_2d(aes(fill = ..count..), geom = "raster", contour = FALSE) +
  scale_fill_gradient2(low = "white", high = "darkblue") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + 
  theme_bw() +
  ggtitle("Subtype mL2/3 - 649 cells") +
  xlab("Cell coverage percent") + ylab("Cell MF") 
ggsave("plots/eda_luo2017_mice/2Ddensityplot_cellMF_vs_cellCov_subtype_mL23_649cells.png", width = 6.5, height = 5)

values(cells_mL23.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.05)))
  ) %>%
  group_by(cell_cov_quantile) %>%
  summarise(percent_cell_MF_less_than_0.5 = sum(cell_MF < 0.5)/length(cell_MF)) %>%
  ggplot(aes(cell_cov_quantile, percent_cell_MF_less_than_0.5)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  theme_bw() + 
  ggtitle("Subtype mL2/3 - 649 cells") +
  xlab("Cell coverage quantile") + ylab("Proportion of cells MF<0.5") + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_luo2017_mice/barplot_propMFless0.5_vs_QcellCov_subtype_mL23_649cells.png", width = 8, height = 5)

values(cells_mL23.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  filter(cell_cov_percent <= 0.15) %>%
  mutate(MF_range = ifelse(cell_MF < 0.5, "< 0.5", ">= 0.5")) %>%
  ggplot(aes(x = cell_cov, color = MF_range, fill = MF_range)) +
  geom_histogram(aes(y = ..density..), position="identity", binwidth = 2, alpha = 0.5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  ggtitle("Subtype mL2/3 - 649 cells") +
  xlab("Cell coverage") + ylab("Density")
ggsave("plots/eda_luo2017_mice/1Ddensityplot_cellCov_groupByMFRange_subtype_mL23_649cells.png", width = 8, height = 5)

with(values(cells_mL23.se), 
     hist(cell_MF[order(cell_cov)[1:100000]], freq = T)
)
with(values(cells_mL23.se), 
     hist(cell_MF[cell_cov>=20], freq = T)
)


# ==== pre-QC mL6-2 686 cells ====
cells_mL62.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_mice/subtype_mL6-2_686cells/"))
values(cells_mL62.se)$cell_cov <- rowSums(assays(cells_mL62.se)$Cov > 0, na.rm = T)
values(cells_mL62.se)$cell_meth <- rowSums(assays(cells_mL62.se)$M/assays(cells_mL62.se)$Cov >= 0.5, na.rm = T) 
values(cells_mL62.se)$cell_MF <- values(cells_mL62.se)$cell_meth/values(cells_mL62.se)$cell_cov

# hist(values(cells.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_mL62.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_mL62.se)$cell_cov, xlim = c(0,100), breaks = 300, freq = F)
quantile(values(cells_mL62.se)$cell_cov)
quantile(values(cells_mL62.se)$cell_cov, probs = 0.05)
values(cells_mL62.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.05)))
  ) %>%
  ggplot(aes(cell_cov_quantile, cell_MF)) +
  geom_boxplot() +
  theme_bw() + 
  ggtitle("Subtype mL6-2 - 686cells") +
  xlab("Cell coverage quantile") + ylab("Cell MF") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_luo2017_mice/boxplot_cellMF_vs_QcellCov_subtype_mL6-2_686cells.png", width = 8, height = 5)

values(cells_mL62.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  slice_sample(prop = 0.1) %>% # suubset a random 10% of the sites to quickly plot
  filter(cell_cov_percent <= 0.15) %>%
  ggplot(aes(cell_cov_percent, cell_MF)) +
  stat_density_2d(aes(fill = ..count..), geom = "raster", contour = FALSE) +
  scale_fill_gradient2(low = "white", high = "darkblue") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + 
  theme_bw() +
  ggtitle("Subtype mL6-2 - 686cells") +
  xlab("Cell coverage percent") + ylab("Cell MF") 
ggsave("plots/eda_luo2017_mice/2Ddensityplot_cellMF_vs_cellCov_subtype_mL6-2_686cells.png", width = 6.5, height = 5)


values(cells_mL62.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.05)))
  ) %>%
  group_by(cell_cov_quantile) %>%
  summarise(percent_cell_MF_less_than_0.5 = sum(cell_MF < 0.5)/length(cell_MF)) %>%
  ggplot(aes(cell_cov_quantile, percent_cell_MF_less_than_0.5)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  theme_bw() + 
  ggtitle("Subtype mL6-2 - 686cells") +
  xlab("Cell coverage quantile") + ylab("Proportion of cells MF<0.5") + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_luo2017_mice/barplot_propMFless0.5_vs_QcellCov_subtype_mL6-2_686cells.png", width = 8, height = 5)

values(cells_mL62.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  filter(cell_cov_percent <= 0.15) %>%
  mutate(MF_range = ifelse(cell_MF < 0.5, "< 0.5", ">= 0.5")) %>%
  ggplot(aes(x = cell_cov, color = MF_range, fill = MF_range)) +
  geom_histogram(aes(y = ..density..), position="identity", binwidth = 2, alpha = 0.5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  ggtitle("Subtype mL6-2 - 686cells") +
  xlab("Cell coverage") + ylab("Density")
ggsave("plots/eda_luo2017_mice/1Ddensityplot_cellCov_groupByMFRange_subtype_mL6-2_686cells.png", width = 8, height = 5)

with(values(cells_mL62.se), 
     hist(cell_MF[order(cell_cov)[1:100000]], freq = T)
)
with(values(cells_mL62.se), 
     hist(cell_MF[cell_cov>=20], freq = T)
)



# ==== pre-QC mPv 136 cells ====
cells_mPv.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_luo2017_mice/subtype_mPv_136cells/"))
values(cells_mPv.se)$cell_cov <- rowSums(assays(cells_mPv.se)$Cov > 0, na.rm = T)
values(cells_mPv.se)$cell_meth <- rowSums(assays(cells_mPv.se)$M/assays(cells_mPv.se)$Cov >= 0.5, na.rm = T) 
values(cells_mPv.se)$cell_MF <- values(cells_mPv.se)$cell_meth/values(cells_mPv.se)$cell_cov

# hist(values(cells.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_mPv.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_mPv.se)$cell_cov, xlim = c(0,100), breaks = 300, freq = F)
quantile(values(cells_mPv.se)$cell_cov)
quantile(values(cells_mPv.se)$cell_cov, probs = 0.05)
values(cells_mPv.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.1)))
  ) %>%
  ggplot(aes(cell_cov_quantile, cell_MF)) +
  geom_boxplot() +
  theme_bw() + 
  ggtitle("Subtype mPv - 136 cells") +
  xlab("Cell coverage quantile") + ylab("Cell MF") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_luo2017_mice/boxplot_cellMF_vs_QcellCov_subtype_mPv_136cells.png", width = 8, height = 5)

values(cells_mPv.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  slice_sample(prop = 0.1) %>% # suubset a random 10% of the sites to quickly plot
  filter(cell_cov_percent <= 0.15) %>%
  ggplot(aes(cell_cov_percent, cell_MF)) +
  stat_density_2d(aes(fill = ..count..), geom = "raster", contour = FALSE) +
  scale_fill_gradient2(low = "white", high = "darkblue") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + 
  theme_bw() +
  ggtitle("Subtype mPv - 136 cells") +
  xlab("Cell coverage percent") + ylab("Cell MF") 
ggsave("plots/eda_luo2017_mice/2Ddensityplot_cellMF_vs_cellCov_subtype_mPv_136cells.png", width = 6.5, height = 5)


values(cells_mPv.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.1)))
  ) %>%
  group_by(cell_cov_quantile) %>%
  summarise(percent_cell_MF_less_than_0.5 = sum(cell_MF < 0.5)/length(cell_MF)) %>%
  ggplot(aes(cell_cov_quantile, percent_cell_MF_less_than_0.5)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  theme_bw() + 
  ggtitle("Subtype mPv - 136 cells") +
  xlab("Cell coverage quantile") + ylab("Proportion of cells MF<0.5") + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_luo2017_mice/barplot_propMFless0.5_vs_QcellCov_subtype_mPv_136cells.png", width = 8, height = 5)

values(cells_mPv.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  filter(cell_cov_percent <= 0.15) %>%
  mutate(MF_range = ifelse(cell_MF < 0.5, "< 0.5", ">= 0.5")) %>%
  ggplot(aes(x = cell_cov, color = MF_range, fill = MF_range)) +
  geom_histogram(aes(y = ..density..), position="identity", binwidth = 1, alpha = 0.5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  ggtitle("Subtype mPv - 136 cells") +
  xlab("Cell coverage") + ylab("Density")
ggsave("plots/eda_luo2017_mice/1Ddensityplot_cellCov_groupByMFRange_subtype_mPv_136cells.png", width = 8, height = 5)

with(values(cells_mPv.se), 
     hist(cell_MF[order(cell_cov)[1:100000]], freq = T)
)
with(values(cells_mPv.se), 
     hist(cell_MF[cell_cov>=20], freq = T)
)




# ==== Overlap of low/high coverage sites across subtypes ==== 

rtail_mL23.gr <- granges(cells_mL23.se)[granges(cells_mL23.se)$cell_cov > quantile(granges(cells_mL23.se)$cell_cov, 0.9)]
ltail_mL23.gr <- granges(cells_mL23.se)[granges(cells_mL23.se)$cell_cov <= quantile(granges(cells_mL23.se)$cell_cov, 0.1)]
rtail_mL62.gr <- granges(cells_mL62.se)[granges(cells_mL62.se)$cell_cov > quantile(granges(cells_mL62.se)$cell_cov, 0.9)]
ltail_mL62.gr <- granges(cells_mL62.se)[granges(cells_mL62.se)$cell_cov <= quantile(granges(cells_mL62.se)$cell_cov, 0.1)]
rtail_mPv.gr <- granges(cells_mPv.se)[granges(cells_mPv.se)$cell_cov > quantile(granges(cells_mPv.se)$cell_cov, 0.9)]
ltail_mPv.gr <- granges(cells_mPv.se)[granges(cells_mPv.se)$cell_cov <= quantile(granges(cells_mPv.se)$cell_cov, 0.1)]

table(countOverlaps(rtail_mL23.gr, rtail_mL62.gr)) / length(rtail_mL23.gr)
table(countOverlaps(rtail_mPv.gr, rtail_mL62.gr)) / length(rtail_mPv.gr)

table(countOverlaps(ltail_mL23.gr, ltail_mL62.gr)) / length(ltail_mL23.gr)
table(countOverlaps(ltail_mPv.gr, ltail_mL62.gr)) / length(ltail_mPv.gr)





















