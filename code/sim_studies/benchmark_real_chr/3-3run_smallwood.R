# ref: Smallwood, S., Lee, H., Angermueller, C. et al. https://doi.org/10.1038/nmeth.3035
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

NV <- 2000


for (N in c(2000)) {
  # for (NP in c(2,3,4,5)) {
  for (NP in c(8,12,20)) {
    for (sparseLevel in 1:3) {
      
      cat("N =", N, "NP =", NP, "\n")
      
      subtype <- "IT-L23_Cux1"
      chromosome <- "chr1"
      seed <- 2022
      alpha <- 0.05
      
      # load input
      dir <- paste0("data/interim/sim_studies/benchmark_real_chr/smallwood/input/pseudoChr_",
                    subtype, "_", chromosome, "_",
                    N, "cells_", NP, "subpops_",
                    NV, "VMRs_sparseLevel", sparseLevel, 
                    "_seed", seed)
      SE <- loadHDF5SummarizedExperiment(dir)
      
      # run model
      t1 <- proc.time()
      
      r_hat_ij <- (assays(SE)$M + 1) / (assays(SE)$Cov + 2)
      var_r_hat_ij <- sqrt(r_hat_ij*(1-r_hat_ij) / assays(SE)$Cov)
      
      w_ij <- 1 / var_r_hat_ij
      w_i <- rowSums(w_ij)
      r_bar_i <- rowSums(r_hat_ij * w_ij) / w_i
      
      n_i_sq <- w_i / (w_i^2 - rowSums(w_ij^2))
      v_hat_i <- n_i_sq * rowSums(w_ij * (r_hat_ij - r_bar_i)^2)
      
      n_i <- sqrt(n_i_sq)
      chi_sq_i <- map_dbl(1:length(n_i), ~ qchisq(p = 1-alpha/2, df = n_i[.x]))
      
      v_hat_i_lb <- n_i * v_hat_i / chi_sq_i
      
      t2 <- proc.time()
      
      # record time elapsed
      time <- round((t2 - t1)[3]/60, 2)
      fwrite(
        data.frame(time = time, method = "smallwood", n_cores = NA),
        paste0("data/interim/sim_studies/benchmark_real_chr/smallwood/output/modelTime_",
               N, "cells_", NP, "subpops_",
               NV, "VMRs_sparseLevel", sparseLevel, 
               "_seed", seed, ".txt")
      )
      
      values(SE)$var_lb <- v_hat_i_lb
      # save model output
      saveRDS(
        granges(SE), 
        paste0(
          "data/interim/sim_studies/benchmark_real_chr/smallwood/output/pseudoChr_",
          subtype, "_", chromosome, "_",
          N, "cells_", NP, "subpops_",
          NV, "VMRs_sparseLevel", sparseLevel, 
          "_seed", seed, "_varLowerBound.rds"
        )
      )
    }
  }
  
  cat("\n\n")
  
}
