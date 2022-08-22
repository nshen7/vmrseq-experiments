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
# for (N in c(500)) {
# for (N in c(1000)) {
  # for (NP in c(2,3,4,5,8,12,20)) {
  for (NP in c(4)) { 
    cat("N =", N, "NP =", NP, "\n")
    
    subtype <- "IT-L23_Cux1"
    chromosome <- "chr1"
    seed <- 2022
    
    # load input
    dir <- paste0("data/interim/sim_studies/benchmark_sim_chr/simulated/simChr_",
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
      paste0("data/interim/sim_studies/benchmark_sim_chr/vmrseq/output/modelTime_", 
             N, "cells_", NP, "subpops_", 
             NV, "VMRs.txt")
    )
    
    # save model output
    saveRDS(res, paste0(
      "data/interim/sim_studies/benchmark_sim_chr/vmrseq/output/simChr_",
      subtype, "_", chromosome, "_", 
      N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed, "_vmrseqOutput.rds"
    ))
    
    cat("\n\n")
  }
}



# my_eval <- function(res_vmrseq, penalty) {
#   res_gr <- res_vmrseq$gr_qced
#   index <- which(res_gr$loglik_diff < penalty) # indices of sites failed to pass penalty threshold
#   res_gr$vmr_index[index] <- NA
#   res_gr$vmr_num_cpg[index] <- NA
#     
#   smr_true <- res_gr %>% data.frame %>% 
#     filter(is_vml) %>%
#     group_by(vmr_name) %>%
#     summarise(n_cpg = n(), 
#               n_dectected = sum(!is.na(vmr_index)), 
#               n_incr = sum(!is.na(cr_index)),
#               cr_name = max(cr_index, na.rm = T),
#               pi1 = max(pi1, na.rm = T),
#               med_total = median(total),
#               loglik_diff = max(loglik_diff, na.rm = T))
#   power <- sum(smr_true$n_dectected > 0) / nrow(smr_true)
#   cat("Power =", power, "; ")
#   smr_detected <- res_gr %>% data.frame %>% 
#     filter(!is.na(vmr_index)) %>%
#     group_by(vmr_index) %>%
#     summarise(cr_index = max(cr_index),
#               n_cpg = n(), 
#               n_true = sum(!is.na(vmr_name)), 
#               pi1 = max(pi1, na.rm = T),
#               optim_pi = max(pi),
#               mean_MF = mean(meth/total),
#               loglik_diff = max(loglik_diff))
#   fdr <- 1-sum(smr_detected$n_true > 0) / nrow(smr_detected)
#   cat("FDR =", fdr, "\n")
#   return(list(power = power, fdr = fdr, true = smr_true, detected = smr_detected))
# }
# 
# for (penalty in 1:10) {
#   smr <- my_eval(res, penalty = penalty)
# }
# 
