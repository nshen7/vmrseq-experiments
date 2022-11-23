.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(scales)

read_dir <- "data/interim/case_studies/luo2017mice_subset_hom/result_summary/"

plot_dir <- "plots/case_studies/luo2017mice_subset_hom/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

# ==== read in results from various methods ====
sites.gr <- readRDS(paste0(read_dir, "cpg_sites.rds"))
res <- list(
  vseq = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs")),
  vseq_cr = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs")),
  scbs = loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs")),
  smwd = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs")),
  scmet = SummarizedExperiment(rowRanges=GRangesList())
)

methods <- c("vmrseq", "CR in vmrseq", "scbs", "smallwood", "scMET")
methods <- factor(methods, levels = methods)

# ==== plots ====
# Basic settings
colors <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = colors[1], "CR in vmrseq" = colors[2], 
                  "scbs" = colors[3], "smallwood" = colors[4], "scMET" = colors[5])

# Number of sites in detected regions
n_sites <- sapply(res, function(se) findOverlaps(granges(se), sites.gr) %>% length())
ggplot() + 
  geom_bar(aes(x = methods, y = n_sites, color = methods, fill = methods), stat = "identity") + 
  scale_y_continuous(labels = scientific, limits = c(0, 2.5e6), name = "Total number of CpGs in detected regions") +
  scale_color_manual(values = COLORVALUES) + 
  scale_fill_manual(values = COLORVALUES) + 
  xlab("Methods") + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_line(color = "light grey"))
ggsave(paste0(plot_dir, "barplot_nSites_vs_methods.png"), height = 4, width = 4)  

# Percentage of sites in detected regions
pct <- n_sites / length(sites.gr)
ggplot() + 
  geom_bar(aes(x = methods, y = pct, color = methods, fill = methods), stat = "identity") + 
  scale_y_continuous(labels = percent, limits = c(0, 0.118), name = "Percentage of CpGs in detected regions") +
  scale_color_manual(values = COLORVALUES) + 
  scale_fill_manual(values = COLORVALUES) + 
  xlab("Methods") + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_line(color = "light grey"))
ggsave(paste0(plot_dir, "barplot_pctSites_vs_methods.png"), height = 4, width = 4)  


