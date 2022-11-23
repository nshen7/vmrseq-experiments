.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
n_cores <- 22
register(MulticoreParam(workers = n_cores))

read_dir <- "data/interim/case_studies/luo2017mice_subset_het/vmrseq/input/"
write_dir <- "data/interim/case_studies/luo2017mice_subset_het/vmrseq/output/"

chrs <- paste0("chr", 1:19)
# for (chr in chrs[1:7]) {
for (chr in chrs[8:19]) {
  # load input
  SE <- loadHDF5SummarizedExperiment(paste0(read_dir, chr))
  
  # QC: remove sites with across-cell coverage < 3
  total <- rowSums(assays(SE)[[1]]>=0, na.rm = T)
  SE <- subset(SE, total >= 3)
  
  # run vmrseq
  gr <- vmrseq.smooth(SE)
  fit <- vmrseq.fit(gr, alpha = 0.05)
  
  # save model output
  saveRDS(fit, paste0(write_dir, chr, ".rds"))
}
