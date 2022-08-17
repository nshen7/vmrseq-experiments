.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
devtools::load_all("../vmrseq-package/vmrseq/")

for (N in c(100)) {
# for (N in c(500)) {
# for (N in c(1000)) {
  for (NP in c(2,3,4,5,8,12,20)) {
    subtype <- "IT-L23_Cux1"
    chromosome <- "chr1"
    seed <- 2022
    
    dir <- paste0("data/interim/sim_studies/benchmark_sim_chr/simulated/simChr_",
                  subtype, "_", chromosome, "_", 
                  N, "cells_", NP, "subpops", "_seed", seed)
    gr <- loadHDF5SummarizedExperiment(dir) %>% granges

    res <- vmrseq(gr, penalty = 0)
    saveRDS(res, paste0(
      "data/interim/sim_studies/benchmark_sim_chr/vmrseq/output/simChr_",
      subtype, "_", chromosome, "_", 
      N, "cells_", NP, "subpops", "_seed", seed, "_vmrseqOutput.rds"
    ))
  }
}