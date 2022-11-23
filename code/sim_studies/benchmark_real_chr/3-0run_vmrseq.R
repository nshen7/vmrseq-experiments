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
  # for (NP in c(12)) {
  # for (NP in c(20)) {
  # for (NP in c(2)) {
  # for (NP in c(3)) {
  # for (NP in c(4)) {
  # for (NP in c(5)) { 
  for (NP in c(8)) {
    # for (sparseLevel in 1) {
    # for (sparseLevel in 2) {
    for (sparseLevel in 3) {
      cat("N =", N, "NP =", NP, "\n")
      cat("alpha = ")

      subtype <- "IT-L23_Cux1"
      chromosome <- "chr1"
      seed <- 2022

      # load input
      dir <- paste0("data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
                    subtype, "_", chromosome, "_",
                    N, "cells_", NP, "subpops_",
                    NV, "VMRs_sparseLevel", sparseLevel,
                    "_seed", seed)
      SE <- loadHDF5SummarizedExperiment(dir)

      # QC: remove sites with across-cell coverage < 3
      total <- rowSums(assays(SE)[[1]]>=0, na.rm = T)
      SE <- subset(SE, total >= 3)

      gr <- vmrseq.smooth(SE)
      for (alpha in c(0.1)){
      # for (alpha in c(seq(0.001, 0.005, 0.001),
      #                 seq(0.01, 0.1, 0.01))) {
      # for (alpha in c(0.12, 0.15, 0.2, 0.3, 0.4)) {
        # run model
        t1 <- proc.time()
        fit <- vmrseq.fit(gr, alpha)
        t2 <- proc.time()

        # record time elapsed
        time <- round((t2 - t1)[3]/60, 2)
        fwrite(
          data.frame(time = time, method = "vmrseq", n_cores = n_cores),
          paste0("data/interim/sim_studies/benchmark_real_chr/vmrseq/output/modelTime_",
                 N, "cells_", NP, "subpops_",
                 NV, "VMRs_sparseLevel", sparseLevel,
                 "_alpha", alpha, ".txt")
        )

        # save model output
        saveRDS(fit, paste0(
          "data/interim/sim_studies/benchmark_real_chr/vmrseq/output/pseudoChr_",
          subtype, "_", chromosome, "_",
          N, "cells_", NP, "subpops_",
          NV, "VMRs_sparseLevel", sparseLevel,
          "_alpha", alpha, "_seed", seed, "_vmrseqOutput.rds"
        ))

        cat(paste0(alpha, ", "))
      }
    }
    cat("\n\n")
  }
}


# for (N in c(200)) {
#   # for (N in c(500, 2000)) {
#   # for (NP in c(2,3,4,5,8,12,20)) {
#   for (NP in c(4)) {
#     # for (sparseLevel in 1:3) {
#       for (sparseLevel in 1) {
#       cat("N =", N, "NP =", NP, "\n")
#       cat("alpha = ")
#       
#       subtype <- "IT-L23_Cux1"
#       chromosome <- "chr1"
#       seed <- 2022
#       
#       # load input
#       dir <- paste0("data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
#                     subtype, "_", chromosome, "_",
#                     N, "cells_", NP, "subpops_",
#                     NV, "VMRs_sparseLevel", sparseLevel, 
#                     "_seed", seed)
#       SE <- loadHDF5SummarizedExperiment(dir)
#       
#       # QC: remove sites with across-cell coverage < 3
#       total <- rowSums(assays(SE)[[1]]>=0, na.rm = T)
#       SE <- subset(SE, total >= 3)
#       
#       gr <- vmrseq.smooth(SE, meanSmooth = TRUE, bpSpan = 500)
#       # for (alpha in c(seq(0.001, 0.005, 0.001), seq(0.01, 0.04, 0.01))) {
#       # for (alpha in c(seq(0.05, 0.1, 0.01))) {
#       # for (alpha in c(0.12, 0.15, 0.2)) {
#       for (alpha in c(0.3, 0.4)) {
#         # run model
#         t1 <- proc.time()
#         fit <- vmrseq.fit(gr, alpha)
#         t2 <- proc.time()
#         
#         # record time elapsed
#         time <- round((t2 - t1)[3]/60, 2)
#         fwrite(
#           data.frame(time = time, method = "vmrseq", n_cores = n_cores),
#           paste0("data/interim/sim_studies/benchmark_real_chr/vmrseq/output/modelTime_",
#                  N, "cells_", NP, "subpops_",
#                  NV, "VMRs_sparseLevel", sparseLevel, 
#                  "_alpha", alpha, "_meanSmoothed.txt")
#         )
#         
#         # save model output
#         saveRDS(fit, paste0(
#           "data/interim/sim_studies/benchmark_real_chr/vmrseq/output/pseudoChr_",
#           subtype, "_", chromosome, "_",
#           N, "cells_", NP, "subpops_",
#           NV, "VMRs_sparseLevel", sparseLevel, 
#           "_alpha", alpha, "_seed", seed, "_vmrseqOutput_meanSmoothed.rds"
#         ))
#         
#         cat(paste0(alpha, ", "))
#       }
#     }
#     cat("\n\n")
#   }
# }

