.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
n_cores <- 30
register(MulticoreParam(workers = n_cores))

NV <- 2000

for (N in c(200)) {
# for (N in c(100, 500)) {
# for (N in c(500)) {
# for (N in c(1000)) {
  # for (NP in c(2,3,4,5,8,12,20)) {
  for (NP in c(4)) {
    cat("N =", N, "NP =", NP, "\n")
    cat("qVar = ")

    subtype <- "IT-L23_Cux1"
    chromosome <- "chr1"
    seed <- 2022

    # load input
    dir <- paste0("data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
                  subtype, "_", chromosome, "_",
                  N, "cells_", NP, "subpops_",
                  NV, "VMRs_seed", seed)
    SE <- loadHDF5SummarizedExperiment(dir)

    # QC: remove sites with across-cell coverage < 3
    total <- rowSums(assays(SE)[[1]]>=0, na.rm = T)
    SE <- subset(SE, total >= 3)

    gr <- vmrseq.smooth(SE)
    # for (qVar in c(0.005, seq(0.01, 0.1, 0.01))) {
    # for (qVar in c(0.12, 0.15, 0.2)) {
    for (qVar in c(0.3, 0.4)) {
      # run model
      t1 <- proc.time()
      fit <- vmrseq.fit(gr, cutoff = quantile(gr$var, prob = 1-qVar))
      t2 <- proc.time()

      # record time elapsed
      time <- round((t2 - t1)[3]/60, 2)
      fwrite(
        data.frame(time = time, method = "vmrseq", n_cores = n_cores),
        paste0("data/interim/sim_studies/benchmark_real_chr/vmrseq/output/modelTime_",
               N, "cells_", NP, "subpops_",
               NV, "VMRs_qVar", qVar, ".txt")
      )

      # save model output
      saveRDS(fit, paste0(
        "data/interim/sim_studies/benchmark_real_chr/vmrseq/output/pseudoChr_",
        subtype, "_", chromosome, "_",
        N, "cells_", NP, "subpops_",
        NV, "VMRs_qVar", qVar, "_seed", seed, "_vmrseqOutput.rds"
      ))

      cat(paste0(qVar, ", "))
    }

    cat("\n\n")
  }
}


# for (N in c(200)) {
#   for (NP in c(4)) {
#     cat("N =", N, "NP =", NP, "\n")
#     cat("qVar = ")
# 
#     subtype <- "IT-L23_Cux1"
#     chromosome <- "chr1"
#     seed <- 2022
# 
#     # load input
#     dir <- paste0("data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
#                   subtype, "_", chromosome, "_",
#                   N, "cells_", NP, "subpops_",
#                   NV, "VMRs_seed", seed)
#     SE <- loadHDF5SummarizedExperiment(dir)
# 
#     # QC: remove sites with across-cell coverage < 3
#     total <- rowSums(assays(SE)[[1]]>=0, na.rm = T)
#     SE <- subset(SE, total >= 3)
# 
#     # for (qVar in c(0.005, seq(0.01, 0.09, 0.01))) {
#     # for (qVar in c(0.1, 0.12, 0.15, 0.2)) {
#     for (qVar in c(0.3, 0.4)) {
#       # run model
#       t1 <- proc.time()
#       res <- vmrseq(SE, qVar = qVar, acrsSmooth = FALSE)
#       t2 <- proc.time()
# 
#       # record time elapsed
#       time <- round((t2 - t1)[3]/60, 2)
#       fwrite(
#         data.frame(time = time, method = "vmrseq", n_cores = n_cores),
#         paste0("data/interim/sim_studies/benchmark_real_chr/vmrseq/output/modelTime_",
#                N, "cells_", NP, "subpops_",
#                NV, "VMRs_qVar", format(qVar, scientific = FALSE), "_noAcrsSmooth.txt")
#       )
# 
#       # save model output
#       saveRDS(res, paste0(
#         "data/interim/sim_studies/benchmark_real_chr/vmrseq/output/pseudoChr_",
#         subtype, "_", chromosome, "_",
#         N, "cells_", NP, "subpops_",
#         NV, "VMRs_qVar", qVar, "_seed", seed, "_vmrseqOutput_noAcrsSmooth.rds"
#       ))
# 
#       cat(paste0(qVar, ", "))
#     }
# 
#     cat("\n\n")
#   }
# }




