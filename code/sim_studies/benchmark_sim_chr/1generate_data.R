.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
source("code/sim_studies/benchmark_sim_chr/helper_functions_generate_chr.R")
devtools::load_all("../vmrseq-package/vmrseq/")
suppressPackageStartupMessages(library(HDF5Array))

NV <- 2000
out_dir <- paste0(getwd(), "/", "data/interim/sim_studies/benchmark_sim_chr/simulated/")

for (N in c(200)) { 
# for (N in c(100)) { 
  # for (NP in c(2,3,4,5,8,12,20)) {
  for (NP in c(4)) {
    # simulateChr(N, NP, NV, sparseLevel = 1, seed = 2022, out_dir)
    # simulateChr(N, NP, NV, sparseLevel = 2, seed = 2022, out_dir)
    simulateChr(N, NP, NV, sparseLevel = 3, seed = 2022, out_dir)
    cat("N =", N, "; NP =", NP, "\n")
  }
}


