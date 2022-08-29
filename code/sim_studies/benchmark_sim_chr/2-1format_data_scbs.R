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
    MF = assays(cell.se)$M_mat*100,
    meth = assays(cell.se)$M_mat,
    unmeth = 1-assays(cell.se)$M_mat
  ) %>% filter(!is.na(MF))
  
  fwrite(cell.df, 
         file = paste0(folder, "/cell_", i, ".cov"),
         row.names = F, col.names = F, sep = "\t", quote = F)
  return(i)
}

# ==== main ====
subtype <- "IT-L23_Cux1"
chromosome <- "chr1"
seed <- 2022
NV <- 3000

# for (N in c(2000)) {
# for (N in c(1000)) {
# for (N in c(500)) {
for (N in c(200)) {
# for (N in c(100)) {
  # for (NP in c(2,3,4,5,8,12,20)) {
  for (NP in c(4)) {
    # Load raw data
    load_dir <- paste0(
      "data/interim/sim_studies/benchmark_sim_chr/simulated/simChr_",
      subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed
    )
    se <- loadHDF5SummarizedExperiment(dir = load_dir)
    
    # format into scbs input
    folder <- paste0(
      "data/interim/sim_studies/benchmark_sim_chr/scbs/input/simChr_",
      subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed
    )
    if (!file.exists(folder)) dir.create(folder)
    
    invisible(bplapply(1:ncol(se), formatCell, se = se, folder = folder))
    cat("N =", N, "; NP =", NP, "\n")
  }
}
