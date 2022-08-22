.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
n_cores <- 14
register(MulticoreParam(workers = n_cores))

NV <- 5000

for (N in c(200)) {
# for (N in c(100, 500)) {
# for (N in c(500)) {
# for (N in c(1000)) {
  # for (NP in c(2,3,4,5,8,12,20)) {
  for (NP in c(4)) { 
    cat("N =", N, "NP =", NP, "\n")
    
    subtype <- "IT-L23_Cux1"
    chromosome <- "chr1"
    seed <- 2022
    
    # load input
    dir <- paste0("data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
                  subtype, "_", chromosome, "_", 
                  N, "cells_", NP, "subpops_", 
                  NV, "VMRs_seed", seed)
    gr <- loadHDF5SummarizedExperiment(dir) %>% granges

    # run model
    t1 <- proc.time()
    res <- vmrseq(gr, penalty = 0)
    t2 <- proc.time()
    
    # record time elapsed
    time <- round((t2 - t1)[3]/60, 2)
    fwrite(
      data.frame(time = time, method = "vmrseq", n_cores = n_cores),  
      paste0("data/interim/sim_studies/benchmark_real_chr/vmrseq/output/modelTime_", 
             N, "cells_", NP, "subpops.txt")
    )
    
    # save model output
    saveRDS(res, paste0(
      "data/interim/sim_studies/benchmark_real_chr/vmrseq/output/pseudoChr_",
      subtype, "_", chromosome, "_", 
      N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed, "_vmrseqOutput.rds"
    ))
    
    for (penalty in 1:20) {
      cat("penalty =", penalty, "\n")
      smr <- my_eval(res, penalty = penalty)
    }
    
    cat("\n\n")
  }
}



