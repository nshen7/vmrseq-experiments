source("code/SETPATHS.R")
source("code/sim_studies/benchmark_sim_chr/helper_functions_generate_chr.R")
devtools::load_all("../vmrseq-package/vmrseq/")
suppressPackageStartupMessages(library(HDF5Array))

NV <- 2000
out_dir <- paste0(getwd(), "/", "data/interim/sim_studies/benchmark_sim_chr/simulated/")

for (N in c(2000)) {
  # for (NP in c(2)) {
  # for (NP in c(3)) {
  # for (NP in c(4)) {
  # for (NP in c(5)) {
  # for (NP in c(8)) {
  for (NP in c(12)) {
  # for (NP in c(20)) {
    # simulateChr(N, NP, NV, sparseLevel = 1, seed = 2022, out_dir)
    # simulateChr(N, NP, NV, sparseLevel = 2, seed = 2022, out_dir)
    simulateChr(N, NP, NV, sparseLevel = 3, seed = 2022, out_dir)
    cat("N =", N, "; NP =", NP, "\n")
  }
}


