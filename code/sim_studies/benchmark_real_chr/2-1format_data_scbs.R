.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") 
library(tidyverse)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)

# ==== util ====
formatCell <- function(i, se, folder) {
  
  cell.se <- se[, i]
  cell.df <- data.frame(
    seqnames = seqnames(cell.se),
    start = start(cell.se),
    end = end(cell.se),
    MF = assays(cell.se)$MF*100,
    meth = assays(cell.se)$MF,
    unmeth = 1-assays(cell.se)$MF
  ) %>% filter(!is.na(MF))
  
  fwrite(cell.df, 
         file = paste0(folder, "/cell_", i, ".cov"),
         row.names = F, col.names = F, sep = "\t", quote = F)
  return(i)
}

# ==== main ====

for (N in c(100, 500, 1000, 2000)) {
  for (NP in c(2,3,4,5,8,12,20)) {
# for (N in c(5000)) {
  # for (NP in c(4,5)) { 
  # for (NP in c(8)) { # TODO
    
    # Load raw data
    load_dir <- paste0(
      "data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
      subtype, "_", chromosome, "_", 
      N, "cells_", NP, "subpops", "_seed", seed
    )
    se <- loadHDF5SummarizedExperiment(dir = load_dir)
    
    # format into scbs input
    folder <- paste0(
      "data/interim/sim_studies/benchmark_real_chr/scbs/pseudoChr_",
      subtype, "_", chromosome, "_", 
      N, "cells_", NP, "subpops", "_seed", seed
    )
    if (!file.exists(folder)) dir.create(folder)
    
    invisible(bplapply(1:ncol(se), formatCell, se = se, folder = folder))
    cat("N =", N, "; NP =", NP, "\n")
  }
}
