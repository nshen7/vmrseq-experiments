source("code/SETPATHS.R")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
library(recommenderlab)
library(Matrix)

# ==== util ====
formatCell <- function(i, se, folder) {
  cell.se <- se[, i]
  M_vec <- dropNA2matrix(as(assays(cell.se)$M_mat, "dgCMatrix")) %>% as.vector()
  cell.df <- data.frame(
    seqnames = seqnames(cell.se),
    start = start(cell.se),
    end = end(cell.se),
    MF = M_vec*100,
    meth = M_vec,
    unmeth = 1-M_vec
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
NV <- 2000

for (N in c(2000)) {
  for (NP in c(8)) {
  # for (NP in c(2,3,4,5,8,12,20)) {
    for (sparseLevel in 3) {
    # for (sparseLevel in 1:3) {
      # Load raw data
      load_dir <- paste0(
        "data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
        subtype, "_", chromosome, "_", 
        N, "cells_", NP, "subpops_", 
        NV, "VMRs_sparseLevel", sparseLevel, 
        "_seed", seed
      )
      se <- loadHDF5SummarizedExperiment(dir = load_dir)
      
      # format into scbs input
      folder <- paste0(
        "data/interim/sim_studies/benchmark_real_chr/scbs/input/pseudoChr_",
        subtype, "_", chromosome, "_", 
        N, "cells_", NP, "subpops_", 
        NV, "VMRs_sparseLevel", sparseLevel, 
        "_seed", seed
      )
      if (!file.exists(folder)) dir.create(folder)
      
      invisible(bplapply(1:ncol(se), formatCell, se = se, folder = folder))
      
      cat("N =", N, "; NP =", NP, "sparse level = ", sparseLevel, "\n")
    }
  }
}
