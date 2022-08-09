.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
source("code/sim_studies/benchmark_sim_chr/helper_functions_generate_chr.R")
devtools::load_all("../vmrseq-package/vmrseq/")
suppressPackageStartupMessages(library(HDF5Array))

out_dir <- paste0(getwd(), "/", "data/interim/sim_studies/benchmark_sim_chr/simulated/")
for (N in c(5000)) {
# for (N in c(100, 500, 1000)) {
  for (NP in c(2,3,4,5,8,12,20)) {
    se <- generatePseudoChr(N, NV = 1000, NP)     
    saveHDF5SummarizedExperiment(
      se,
      paste0("data/interim/sim_studies/benchmark_sim_chr/simulated/pseudo_chr_",N,"cells_",NP,"subpops/"),
      replace = TRUE
    )
    cat("N =", N, "; NP =", NP, "\n")
  }
}


