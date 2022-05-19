library(data.table)
library(readxl)
library(tidyverse)
library(SummarizedExperiment)
library(GenomicRanges)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
library(bsseq)
library(pheatmap)
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")

# ==== read in metadata ====
metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv")
metadata %>%
  ggplot(aes(Sample, FinalReads/1e6)) +
  geom_boxplot() +
  theme_bw() + ylim(0,5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### high coverage samples
metadata[Sample=="1A_180226", .(.N), by = .(CellClass, SubType)] %>% arrange(desc(N)) %>% head()
metadata[Sample=="8J_190716", .(.N), by = .(CellClass, SubType)] %>% arrange(desc(N)) %>% head()
metadata[Sample=="11F_190214", .(.N), by = .(CellClass, SubType)] %>% arrange(desc(N)) %>% head()
metadata[Sample=="1A_180226", summary(FinalReads)/1e6]
metadata[Sample=="8J_190716", summary(FinalReads)/1e6]
metadata[Sample=="11F_190214", summary(FinalReads)/1e6]

### low coverage samples
metadata[Sample=="2E_180220", .(.N), by = .(CellClass, SubType)] %>% arrange(desc(N)) %>% head()
metadata[Sample=="4A_180205", .(.N), by = .(CellClass, SubType)] %>% arrange(desc(N)) %>% head()
metadata[Sample=="1C_180212", .(.N), by = .(CellClass, SubType)] %>% arrange(desc(N)) %>% head()
metadata[Sample=="2E_180220", summary(FinalReads)/1e6]
metadata[Sample=="4A_180205", summary(FinalReads)/1e6]
metadata[Sample=="1C_180212", summary(FinalReads)/1e6]

# ==== (low coverage 1) pre-QC sample 2E_180220 OLF-Exc Bmpr1b n=564 ====
cells_low1.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample2E_180220_GSE131904_subtype_OLF-Exc_Bmpr1b_564cells/"))
values(cells_low1.se)$cell_cov <- rowSums(assays(cells_low1.se)$Cov > 0, na.rm = T)
values(cells_low1.se)$cell_meth <- rowSums(assays(cells_low1.se)$M/assays(cells_low1.se)$Cov >= 0.5, na.rm = T) 
values(cells_low1.se)$cell_MF <- values(cells_low1.se)$cell_meth/values(cells_low1.se)$cell_cov

# hist(values(cells.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_low1.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_low1.se)$cell_cov, xlim = c(0,100), breaks = 300, freq = F)
quantile(values(cells_low1.se)$cell_cov)
quantile(values(cells_low1.se)$cell_cov, probs = 0.05)
values(cells_low1.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.05)))
         ) %>%
  ggplot(aes(cell_cov_quantile, cell_MF)) +
  geom_boxplot() +
  theme_bw() + 
  ggtitle("Low-coverage sample 2E_180220 (subtype OLF-Exc Bmpr1b - 564 cells)") +
  xlab("Cell coverage quantile") + ylab("Cell MF") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/boxplot_cellMF_vs_QcellCov_lowCov_sample2E_180220_GSE131904_subtype_OLF-Exc_Bmpr1b_564cells.png", width = 8, height = 5)

values(cells_low1.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  slice_sample(prop = 0.1) %>% # subset a random 10% of the sites to quickly plot
  filter(cell_cov_percent <= 0.15) %>%
  ggplot(aes(cell_cov_percent, cell_MF)) +
  stat_density_2d(aes(fill = ..count..), geom = "raster", contour = FALSE) +
  scale_fill_gradient2(low = "white", high = "darkblue") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + 
  theme_bw() +
  ggtitle("Low-coverage sample 2E_180220 (subtype OLF-Exc Bmpr1b - 564 cells)") +
  xlab("Cell coverage percent") + ylab("Cell MF") 
ggsave("plots/eda_liu2021/2Ddensityplot_cellMF_vs_cellCov_lowCov_sample2E_180220_GSE131904_subtype_OLF-Exc_Bmpr1b_564cells.png", width = 6.5, height = 5)

values(cells_low1.se) %>%
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
  ggtitle("Low-coverage sample 2E_180220 (subtype OLF-Exc Bmpr1b - 564 cells)") +
  xlab("Cell coverage quantile") + ylab("Proportion of cells MF<0.5") + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/barplot_propMFless0.5_vs_QcellCov_lowCov_sample2E_180220_GSE131904_subtype_OLF-Exc_Bmpr1b_564cells.png", width = 8, height = 5)

values(cells_low1.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  filter(cell_cov_percent <= 0.15) %>%
  mutate(MF_range = ifelse(cell_MF < 0.5, "< 0.5", ">= 0.5")) %>%
  ggplot(aes(x = cell_cov, color = MF_range, fill = MF_range)) +
  geom_histogram(aes(y = ..density..), position="identity", binwidth = 2, alpha = 0.5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  ggtitle("Low-coverage sample 2E_180220 (subtype OLF-Exc Bmpr1b - 564 cells)") +
  xlab("Cell coverage") + ylab("Density")
ggsave("plots/eda_liu2021/1Ddensityplot_cellCov_groupByMFRange_lowCov_sample2E_180220_GSE131904_subtype_OLF-Exc_Bmpr1b_564cells.png", width = 6.5, height = 5)


granges(cells_low1.se) %>%
  as.data.frame() %>%
  filter(cell_cov <= 0.15 * max(cell_cov)) %>%
  filter(seqnames %in% c("chr1", "chr2")) %>%
  slice_sample(prop = 0.001) %>% # subset a random 10% of the sites to quickly plot
  ggplot(aes(x = start, y = cell_cov, color = cell_MF)) +
  geom_point(alpha = 0.5, size = 0.2) +
  facet_wrap(~seqnames) +
  theme_bw() + 
  viridis::scale_color_viridis(name = "MF", option = "D") +
  ggtitle("Low-coverage sample 2E_180220 (subtype OLF-Exc Bmpr1b - 564 cells)") +
  xlab("Genomic coordinate") + ylab("Cell coverage")
ggsave("plots/eda_liu2021/scatterplot_cellCov_vs_genomCoord_lowCov_sample2E_180220_GSE131904_subtype_OLF-Exc_Bmpr1b_564cells.png", width = 6.5, height = 5)

with(values(cells_low1.se), 
     hist(cell_MF[order(cell_cov)[1:100000]], freq = T)
)
with(values(cells_low1.se), 
     hist(cell_MF[cell_cov>=20], freq = T)
)

# ==== (low coverage 2) pre-QC sample: 4A_180205 Exc IT-L23 Tenm2 190 ====
cells_low2.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample4A_180205_GSE131766_subtype_IT-L23_Tenm2_190cells/"))
values(cells_low2.se)$cell_cov <- rowSums(assays(cells_low2.se)$Cov > 0, na.rm = T)
values(cells_low2.se)$cell_meth <- rowSums(assays(cells_low2.se)$M/assays(cells_low2.se)$Cov >= 0.5, na.rm = T) 
values(cells_low2.se)$cell_MF <- values(cells_low2.se)$cell_meth/values(cells_low2.se)$cell_cov

hist(values(cells_low2.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_low2.se)$cell_cov, breaks = 200, freq = F)
quantile(values(cells_low2.se)$cell_cov)
quantile(values(cells_low2.se)$cell_cov, probs = 0.05)
values(cells_low2.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.1)))
  ) %>%
  ggplot(aes(cell_cov_quantile, cell_MF)) +
  geom_boxplot() +
  theme_bw() + 
  xlab("Cell coverage quantile") + ylab("Cell MF") +
  ggtitle("Low-coverage sample 4A_180205 (subtype IT-L23 Tenm2 Cux1 - 190 cells)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/boxplot_cellMF_vs_QcellCov_lowCov_sample4A_180205_GSE131766_subtype_IT-L23_Tenm2_190cells.png", width = 8, height = 5)

values(cells_low2.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  slice_sample(prop = 0.1) %>% # suubset a random 10% of the sites to quickly plot
  filter(cell_cov_percent <= 0.15) %>%
  ggplot(aes(cell_cov_percent, cell_MF)) +
  stat_density_2d(aes(fill = ..count..), geom = "raster", contour = FALSE) +
  scale_fill_gradient2(low = "white", high = "darkblue") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + 
  ggtitle("Low-coverage sample 4A_180205 (subtype IT-L23 Tenm2 Cux1 - 190 cells)") +
  xlab("Cell coverage percent") + ylab("Cell MF") 
ggsave("plots/eda_liu2021/2Ddensityplot_cellMF_vs_cellCov_lowCov_sample4A_180205_GSE131766_subtype_IT-L23_Tenm2_190cells.png", width = 6.5, height = 5)

values(cells_low2.se) %>%
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
  ggtitle("Low-coverage sample 4A_180205 (subtype IT-L23 Tenm2 Cux1 - 190 cells)") +
  xlab("Cell coverage quantile") + ylab("Proportion of cells MF<0.5") + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/barplot_propMFless0.5_vs_QcellCov_lowCov_sample4A_180205_GSE131766_subtype_IT-L23_Tenm2_190cells.png", width = 8, height = 5)


values(cells_low2.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  filter(cell_cov_percent <= 0.15) %>%
  mutate(MF_range = ifelse(cell_MF < 0.5, "< 0.5", ">= 0.5")) %>%
  ggplot(aes(x = cell_cov, color = MF_range, fill = MF_range)) +
  geom_histogram(aes(y = ..density..), position="identity", binwidth = 1, alpha = 0.5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  ggtitle("Low-coverage sample 4A_180205 (subtype IT-L23 Tenm2 Cux1 - 190 cells)") +
  xlab("Cell coverage") + ylab("Density")
ggsave("plots/eda_liu2021/1Ddensityplot_cellCov_groupByMFRange_lowCov_sample4A_180205_GSE131766_subtype_IT-L23_Tenm2_190cells.png", width = 6.5, height = 5)

with(values(cells_low2.se), 
     hist(cell_MF[order(cell_cov)[1:100000]], freq = T)
)
with(values(cells_low2.se), 
     hist(cell_MF[cell_cov>=quantile(cell_cov, probs = 0.8)], freq = T)
)


# ==== (low-mid coverage 3) pre-QC sample: 1C_180212 Inh OLF Trpc4 n=264 ====
cells_low3.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample1C_180212_GSE131867_subtype_OLF_Trpc4_264cells/"))
values(cells_low3.se)$cell_cov <- rowSums(assays(cells_low3.se)$Cov > 0, na.rm = T)
values(cells_low3.se)$cell_meth <- rowSums(assays(cells_low3.se)$M/assays(cells_low3.se)$Cov >= 0.5, na.rm = T) 
values(cells_low3.se)$cell_MF <- values(cells_low3.se)$cell_meth/values(cells_low3.se)$cell_cov

# hist(values(cells.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_low3.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_low3.se)$cell_cov, xlim = c(0,100), breaks = 300, freq = F)
quantile(values(cells_low3.se)$cell_cov, probs = seq(0,1,0.05))
values(cells_low3.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.1)))
  ) %>%
  ggplot(aes(cell_cov_quantile, cell_MF)) +
  geom_boxplot() +
  theme_bw() + 
  ggtitle("Low-coverage sample 1C_180212 (subtype OLF Trpc4 - 264 cells)") +
  xlab("Cell coverage quantile") + ylab("Cell MF") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/boxplot_cellMF_vs_QcellCov_lowCov_sample1C_180212_GSE131867_subtype_OLF_Trpc4_264cells.png", width = 8, height = 5)

values(cells_low3.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  slice_sample(prop = 0.1) %>% # suubset a random 10% of the sites to quickly plot
  filter(cell_cov_percent <= 0.15) %>%
  ggplot(aes(cell_cov_percent, cell_MF)) +
  stat_density_2d(aes(fill = ..count..), geom = "raster", contour = FALSE) +
  scale_fill_gradient2(low = "white", high = "darkblue") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + 
  ggtitle("Low-coverage sample 1C_180212 (subtype OLF Trpc4 - 264 cells)") +
  xlab("Cell coverage percent") + ylab("Cell MF") 
ggsave("plots/eda_liu2021/2Ddensityplot_cellMF_vs_cellCov_lowCov_sample1C_180212_GSE131867_subtype_OLF_Trpc4_264cellss.png", width = 6.5, height = 5)


values(cells_low3.se) %>%
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
  ggtitle("Low-coverage sample 1C_180212 (subtype OLF Trpc4 - 264 cells)") +
  xlab("Cell coverage quantile") + ylab("Proportion of cells MF<0.5") + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/barplot_propMFless0.5_vs_QcellCov_lowCov_sample1C_180212_GSE131867_subtype_OLF_Trpc4_264cells.png", width = 8, height = 5)

values(cells_low3.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  filter(cell_cov_percent <= 0.15) %>%
  mutate(MF_range = ifelse(cell_MF < 0.5, "< 0.5", ">= 0.5")) %>%
  ggplot(aes(x = cell_cov, color = MF_range, fill = MF_range)) +
  geom_histogram(aes(y = ..density..), position="identity", binwidth = 1, alpha = 0.5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  ggtitle("Low-coverage sample 1C_180212 (subtype OLF Trpc4 - 264 cells)") +
  xlab("Cell coverage") + ylab("Density")
ggsave("plots/eda_liu2021/1Ddensityplot_cellCov_groupByMFRange_lowCov_sample1C_180212_GSE131867_subtype_OLF_Trpc4_264cells.png", width = 6.5, height = 5)

with(values(cells_low3.se), 
     hist(cell_MF[order(cell_cov)[1:100000]], freq = T)
)
with(values(cells_low3.se), 
     hist(cell_MF[cell_cov>=20], freq = T)
)



# ==== (high coverage 1) pre-QC sample: 1A_180226 IT-L23 Cux1 GSE130553 n=506 ====
cells_high1.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample1A_180226_GSE130553_subtype_IT-L23_Cux1_506cells/"))
values(cells_high1.se)$cell_cov <- rowSums(assays(cells_high1.se)$Cov > 0, na.rm = T)
values(cells_high1.se)$cell_meth <- rowSums(assays(cells_high1.se)$M/assays(cells_high1.se)$Cov >= 0.5, na.rm = T) 
values(cells_high1.se)$cell_MF <- values(cells_high1.se)$cell_meth/values(cells_high1.se)$cell_cov

hist(values(cells_high1.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_high1.se)$cell_cov, breaks = 200, freq = F)
quantile(values(cells_high1.se)$cell_cov)
quantile(values(cells_high1.se)$cell_cov, probs = 0.025)
values(cells_high1.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.05)))
  ) %>%
  ggplot(aes(cell_cov_quantile, cell_MF)) +
  geom_boxplot() +
  theme_bw() + 
  xlab("Cell coverage quantile") + ylab("Cell MF") +
  ggtitle("High-coverage sample 1A_180226 (subtype IT-L23 Cux1 - 506 cells)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/boxplot_cellMF_vs_QcellCov_highCov_sample1A_180226_GSE130553_subtype_IT-L23_Cux1_506cells.png", width = 8, height = 5)

values(cells_high1.se) %>%
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
  ggtitle("High-coverage sample 1A_180226 (subtype IT-L23 Cux1 - 506 cells)") +
  xlab("Cell coverage percent") + ylab("Cell MF") 
ggsave("plots/eda_liu2021/2Ddensityplot_cellMF_vs_cellCov_highCov_sample1A_180226_GSE130553_subtype_IT-L23_Cux1_506cells.png", width = 6.5, height = 5)

values(cells_high1.se) %>%
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
  ggtitle("High-coverage sample 1A_180226 (subtype IT-L23 Cux1 - 506 cells)") +
  xlab("Cell coverage quantile") + ylab("Proportion of cells MF<0.5") + ylim(0, 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/barplot_propMFless0.5_vs_QcellCov_highCov_sample1A_180226_GSE130553_subtype_IT-L23_Cux1_506cells.png", width = 8, height = 5)


values(cells_high1.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  filter(cell_cov_percent <= 0.15) %>%
  mutate(MF_range = ifelse(cell_MF < 0.5, "< 0.5", ">= 0.5")) %>%
  ggplot(aes(x = cell_cov, color = MF_range, fill = MF_range)) +
  geom_histogram(aes(y = ..density..), position="identity", binwidth = 2, alpha = 0.5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  ggtitle("High-coverage sample 1A_180226 (subtype IT-L23 Cux1 - 506 cells)") +
  xlab("Cell coverage") + ylab("Density")
ggsave("plots/eda_liu2021/1Ddensityplot_cellCov_groupByMFRange_highCov_sample1A_180226_GSE130553_subtype_IT-L23_Cux1_506cells.png", width = 6.5, height = 5)

with(values(cells_high1.se), 
     hist(cell_MF[cell_cov<=9], freq = T)
)
with(values(cells_high1.se), 
     hist(cell_MF[cell_cov>=quantile(cell_cov, probs = 0.8)], freq = T)
)


# ==== subsampling a high-coverage subtype 1A_180226 to a pseudo-low-coverage sample ====
# cells_high1.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample1A_180226_GSE130553_subtype_IT-L23_Cux1_506cells/"))
# Cov_sub <- as(assays(cells_high1.se)$Cov, "sparseMatrix")
# M_sub <- as(assays(cells_high1.se)$M, "sparseMatrix")
# for(i in 1:ncol(Cov_sub)) {
#   ind <- which(Cov_sub[,i]>0)
#   set.seed(i); ind_rm <- sample(ind, round(length(ind)/2))
#   Cov_sub[ind_rm,i] <- 0
#   M_sub[ind_rm,i] <- 0
#   cat(i, " ")
# }
df <-  data.frame(cell_cov_sub = sapply(values(cells_high1.se)$cell_cov, function(x) rbinom(1, x, 0.5)),
                  cell_MF = values(cells_high1.se)$cell_MF) 
df$cell_meth_sub <- apply(df, 1, function(x) rbinom(1,x[1],x[2]))

df <- df %>%
  as.data.frame() %>%
  filter(cell_cov_sub <= 100 & cell_cov_sub > 0) %>%
  mutate(cell_MF_sub = cell_meth_sub / cell_cov_sub) %>%
  mutate(cell_cov_sub_quantile = cut(cell_cov_sub,
                                     breaks = quantile(c(0, cell_cov_sub), probs = seq(0,1,0.1)))
  ) 

df %>%
  ggplot(aes(cell_cov_sub_quantile, cell_MF)) +
  geom_boxplot() +
  theme_bw() + 
  xlab("Cell coverage quantile") + ylab("Cell MF") +
  ggtitle("Subsampled high-coverage sample 1A_180226 (subtype IT-L23 Cux1 - 506 cells)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





# ==== (high coverage 2) pre-QC sample: 8J_190716  Exc  DG dg-all  GSE158154  n=591 ====
cells_high2.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample8J_190716_GSE158154_subtype_DG_dg-all_591cells/"))
values(cells_high2.se)$cell_cov <- rowSums(assays(cells_high2.se)$Cov > 0, na.rm = T)
values(cells_high2.se)$cell_meth <- rowSums(assays(cells_high2.se)$M/assays(cells_high2.se)$Cov >= 0.5, na.rm = T) 
values(cells_high2.se)$cell_MF <- values(cells_high2.se)$cell_meth/values(cells_high2.se)$cell_cov

hist(values(cells_high2.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_high2.se)$cell_cov, breaks = 200, freq = F)
quantile(values(cells_high2.se)$cell_cov)
quantile(values(cells_high2.se)$cell_cov, probs = 0.05)
values(cells_high2.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.05)))
  ) %>%
  ggplot(aes(cell_cov_quantile, cell_MF)) +
  geom_boxplot() +
  theme_bw() + 
  xlab("Cell coverage quantile") + ylab("Cell MF") +
  ggtitle("High-coverage sample 8J_190716 (subtype DG dg-all - 591 cells)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/boxplot_cellMF_vs_QcellCov_highCov_sample8J_190716_GSE158154_subtype_DG_dg-all_591cells.png", width = 8, height = 5)

values(cells_high2.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  slice_sample(prop = 0.1) %>% # suubset a random 10% of the sites to quickly plot
  filter(cell_cov_percent <= 0.15) %>%
  ggplot(aes(cell_cov_percent, cell_MF)) +
  stat_density_2d(aes(fill = ..count..), geom = "raster", contour = FALSE) +
  scale_fill_gradient2(low = "white", high = "darkblue") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + 
  ggtitle("High-coverage sample 8J_190716 (subtype DG dg-all - 591 cells)") +
  xlab("Cell coverage percent") + ylab("Cell MF") 
ggsave("plots/eda_liu2021/2Ddensityplot_cellMF_vs_cellCov_highCov_sample8J_190716_GSE158154_subtype_DG_dg-all_591cells.png", width = 6.5, height = 5)

values(cells_high2.se) %>%
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
  ggtitle("High-coverage sample 8J_190716 (subtype DG dg-all Cux1 - 591 cells)") +
  xlab("Cell coverage quantile") + ylab("Proportion of cells MF<0.5") + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/barplot_propMFless0.5_vs_QcellCov_highCov_sample8J_190716_GSE158154_subtype_DG_dg-all_591cells.png", width = 8, height = 5)

values(cells_high2.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  filter(cell_cov_percent <= 0.15) %>%
  mutate(MF_range = ifelse(cell_MF < 0.5, "< 0.5", ">= 0.5")) %>%
  ggplot(aes(x = cell_cov, color = MF_range, fill = MF_range)) +
  geom_histogram(aes(y = ..density..), position="identity", binwidth = 2, alpha = 0.5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  ggtitle("High-coverage sample 8J_190716 (subtype DG dg-all Cux1 - 591 cells)") +
  xlab("Cell coverage") + ylab("Density")
ggsave("plots/eda_liu2021/1Ddensityplot_cellCov_groupByMFRange_highCov_sample8J_190716_GSE158154_subtype_DG_dg-all_591cells.png", width = 6.5, height = 5)

with(values(cells_high2.se), 
     hist(cell_MF[order(cell_cov)[1:100000]], freq = T)
)
with(values(cells_high2.se), 
     hist(cell_MF[cell_cov>=quantile(cell_cov, probs = 0.8)], freq = T)
)



# ==== (high coverage 3) pre-QC sample: 11F_190214  Exc  DG dg-all GSE135543 n=871 ====
cells_high3.se <- loadHDF5SummarizedExperiment(dir = paste0("data/processed/processed_liu2021/sample11F_190214_GSE135543_subtype_DG_dg-all_871cells/"))
values(cells_high3.se)$cell_cov <- rowSums(assays(cells_high3.se)$Cov > 0, na.rm = T)
values(cells_high3.se)$cell_meth <- rowSums(assays(cells_high3.se)$M/assays(cells_high3.se)$Cov >= 0.5, na.rm = T) 
values(cells_high3.se)$cell_MF <- values(cells_high3.se)$cell_meth/values(cells_high3.se)$cell_cov

hist(values(cells_high3.se)$cell_MF, breaks = 200, freq = F)
hist(values(cells_high3.se)$cell_cov, breaks = 200, freq = F)
quantile(values(cells_high3.se)$cell_cov)
quantile(values(cells_high3.se)$cell_cov, probs = 0.05)
values(cells_high3.se) %>%
  as.data.frame() %>%
  # filter(cell_cov <= 100) %>%
  mutate(cell_cov_quantile = cut(cell_cov, 
                                 breaks = quantile(c(0,cell_cov), probs = seq(0,1,0.05)))
  ) %>%
  ggplot(aes(cell_cov_quantile, cell_MF)) +
  geom_boxplot() +
  theme_bw() + 
  xlab("Cell coverage quantile") + ylab("Cell MF") +
  ggtitle("High-coverage sample 11F_190214 (subtype DG dg-all - 871 cells)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/boxplot_cellMF_vs_QcellCov_highCov_sample11F_190214_GSE135543_subtype_DG_dg-all_871cells.png", width = 8, height = 5)

values(cells_high3.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  slice_sample(prop = 0.1) %>% # suubset a random 10% of the sites to quickly plot
  filter(cell_cov_percent <= 0.15) %>%
  ggplot(aes(cell_cov_percent, cell_MF)) +
  stat_density_2d(aes(fill = ..count..), geom = "raster", contour = FALSE) +
  scale_fill_gradient2(low = "white", high = "darkblue") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) + 
  ggtitle("High-coverage sample 11F_190214 (subtype DG dg-all - 871 cells)") +
  xlab("Cell coverage percent") + ylab("Cell MF") 
ggsave("plots/eda_liu2021/2Ddensityplot_cellMF_vs_cellCov_highCov_sample11F_190214_GSE135543_subtype_DG_dg-all_871cells.png", width = 6.5, height = 5)


values(cells_high3.se) %>%
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
  ggtitle("High-coverage sample 11F_190214 (subtype DG dg-all - 871 cells)") +
  xlab("Cell coverage quantile") + ylab("Proportion of cells MF<0.5") + ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/eda_liu2021/barplot_propMFless0.5_vs_QcellCov_highCov_sample11F_190214_GSE135543_subtype_DG_dg-all_871cells.png", width = 8, height = 5)

values(cells_high3.se) %>%
  as.data.frame() %>%
  mutate(cell_cov_percent = cell_cov / max(cell_cov)) %>%
  filter(cell_cov_percent <= 0.15) %>%
  mutate(MF_range = ifelse(cell_MF < 0.5, "< 0.5", ">= 0.5")) %>%
  ggplot(aes(x = cell_cov, color = MF_range, fill = MF_range)) +
  geom_histogram(aes(y = ..density..), position="identity", binwidth = 2, alpha = 0.5) +
  theme_bw() + 
  scale_color_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), name = "MF Range") +
  ggtitle("High-coverage sample 11F_190214 (subtype DG dg-all - 871 cells)") +
  xlab("Cell coverage") + ylab("Density")
ggsave("plots/eda_liu2021/1Ddensityplot_cellCov_groupByMFRange_highCov_sample11F_190214_GSE135543_subtype_DG_dg-all_871cells.png", width = 6.5, height = 5)


with(values(cells_high3.se), 
     hist(cell_MF[order(cell_cov)[1:100000]], freq = T)
)
with(values(cells_high3.se), 
     hist(cell_MF[cell_cov>=quantile(cell_cov, probs = 0.8)], freq = T)
)

# ==== total coverage per cell vs. MF ====
total_cov_low1 <- colSums(assays(cells_low1.se)$Cov > 0, na.rm = T)
total_meth_low1 <- colSums(assays(cells_low1.se)$M/assays(cells_low1.se)$Cov >= 0.5, na.rm = T) 
total_cov_low2 <- colSums(assays(cells_low2.se)$Cov > 0, na.rm = T)
total_meth_low2 <- colSums(assays(cells_low2.se)$M/assays(cells_low2.se)$Cov >= 0.5, na.rm = T) 
total_cov_high1 <- colSums(assays(cells_high1.se)$Cov > 0, na.rm = T)
total_meth_high1 <- colSums(assays(cells_high1.se)$M/assays(cells_high1.se)$Cov >= 0.5, na.rm = T) 
total_cov_high2 <- colSums(assays(cells_high2.se)$Cov > 0, na.rm = T)
total_meth_high2 <- colSums(assays(cells_high2.se)$M/assays(cells_high2.se)$Cov >= 0.5, na.rm = T) 

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
data.frame(total_cov = c(total_cov_high1, total_cov_high2, total_cov_low1, total_cov_low2),
           total_meth = c(total_meth_high1, total_meth_high2, total_meth_low1, total_meth_low2),
           cov_level = c(rep("high 1", length(total_cov_high1)), 
                         rep("high 2", length(total_cov_high2)),
                         rep("low 1", length(total_cov_low1)),
                         rep("low 2", length(total_cov_low2))
                         )) %>%
  mutate(total_MF = total_meth / total_cov) %>%
  ggplot(aes(total_cov, total_MF, color = cov_level)) +
  geom_point(size = 1) + 
  geom_smooth(aes(group = cov_level), method = "lm") + 
  scale_color_manual(values = cbp1, name="Sample") +
  theme_bw() +
  xlab("Total covered sites in cell") +
  ylab("mCG/CG") + ylim(0.6,0.9)
ggsave("plots/eda_liu2021/point_MFofSites_vs_totalCoveredSites_4Samples.png", width = 7, height = 5)  


metadata %>%
  filter(`Pass QC` == TRUE) %>%
  ggplot(aes(FinalReads, CG_Frac, color = CellClass)) +
  geom_point(size = 0.5, alpha = 0.5) + 
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("Final reads") +
  ylab("mCG/CG") + ylim(0.6,0.9)
ggsave("plots/eda_liu2021/point_MFofSites_vs_totalCoveredSites_allCells.png", width = 7, height = 5)  

metadata %>%
  filter(`Pass QC` == TRUE) %>%
  ggplot(aes(MappingRate, FinalReads)) +
  geom_point(size = 0.5, alpha = 0.5) + 
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("Mapping Rate") +
  ylab("Final Reads") 

# ==== overlap of low/high coverage sites between subtypes ==== 

### NOTICE: luo data loaded from "code/luo2017_mice_data_processing/cov_vs_MF_eda.R"
switchSubtype <- function(subtype) switch(as.character(subtype),
                                          liu_high1 = granges(cells_high1.se),
                                          liu_high2 = granges(cells_high2.se),
                                          liu_high3 = granges(cells_high3.se),
                                          liu_low1 = granges(cells_low1.se),
                                          liu_low2 = granges(cells_low2.se),
                                          liu_low3 = granges(cells_low3.se),
                                          luo_mL23 = granges(cells_mL23.se),
                                          luo_mL62 = granges(cells_mL62.se),
                                          luo_mPv = granges(cells_mPv.se))

subtypes <- c("luo_mL23", "luo_mL62", "luo_mPv", 
              "liu_high1", "liu_high2", "liu_high3",
              "liu_low1", "liu_low2", "liu_low3")
df <- expand.grid(quest = subtypes, subject = subtypes, stringsAsFactors = F)


ltailExtract <- function(gr_obj) gr_obj[gr_obj$cell_cov <= quantile(gr_obj$cell_cov, 0.1)]
rtailExtract <- function(gr_obj) gr_obj[gr_obj$cell_cov > quantile(gr_obj$cell_cov, 0.9)]

ltailOverlap <- function(quest_subtype, object_subtype) {
  quest <- ltailExtract(switchSubtype(quest_subtype))
  object <- ltailExtract(switchSubtype(object_subtype))
  return(sum(countOverlaps(quest, object)>0) / length(quest))
}
rtailOverlap <- function(quest_subtype, object_subtype) {
  quest <- rtailExtract(switchSubtype(quest_subtype))
  object <- rtailExtract(switchSubtype(object_subtype))
  return(sum(countOverlaps(quest, object)>0) / length(quest))
}

df$ltail_olap_pct <- sapply(1:nrow(df), function(i) ltailOverlap(df$quest[i], df$subject[i]))
df$rtail_olap_pct <- sapply(1:nrow(df), function(i) rtailOverlap(df$quest[i], df$subject[i]))

df %>%
  ggplot(aes(subject, quest, fill = ltail_olap_pct)) + 
  geom_tile() +
  geom_text(aes(label = round(ltail_olap_pct, 2)), size = 3, color = "white") +
  viridis::scale_fill_viridis(limits = c(0,1), name = "Overlap\nPercentage", option = "D") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Overlap of low-coverage sites between subtypes") +
  xlab("subject (reference)")
ggsave("plots/eda_liu2021/heatmap_betweenSubtypes_lowCovSites_olapPct.png", height = 5, width = 6.5)


df %>%
  ggplot(aes(subject, quest, fill = rtail_olap_pct)) + 
  geom_tile() +
  geom_text(aes(label = round(rtail_olap_pct, 2)), size = 3, color = "white") +
  viridis::scale_fill_viridis(limits = c(0,1), name = "Overlap\nPercentage", option = "D") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Overlap of high-coverage sites between subtypes") +
  xlab("subject (reference)")
ggsave("plots/eda_liu2021/heatmap_betweenSubtypes_highCovSites_olapPct.png", height = 5, width = 6.5)



















