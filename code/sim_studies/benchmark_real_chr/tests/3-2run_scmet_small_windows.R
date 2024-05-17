source("code/SETPATHS.R")
library(scMET)
library(tidyverse)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
n_cores <- 22
NV <- 2000
bp_size <- 10000

for (N in c(500)) {
  for (NP in c(8)) {
    for (sparseLevel in 1) {
      
      subtype <- "IT-L23_Cux1"
      chromosome <- "chr1"
      seed <- 2022
      
      # Load input data
      dir <- paste0("data/interim/sim_studies/benchmark_real_chr/scmet/input/pseudoChr_", 
                    subtype, "_", chromosome, "_",
                    bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
                    NV, "VMRs_sparseLevel", sparseLevel, 
                    "_seed", seed, ".txt.gz")
      dt <- fread(dir) %>% 
        mutate(met_reads = as.integer(round(met_reads))) %>%
        dplyr::filter(total_reads >= 3) %>%  # exclude regions that have less than 3 CpGs covered
        group_by(Feature) %>%
        dplyr::filter(n() >= 5) # exclude features that did not have CpG coverage in at least 5 cells
      
      # Run model
      t1 <- proc.time()
      fit_obj <- scmet(Y = dt, L = 4, iter = 50000, seed = 12, n_cores = n_cores)
      t2 <- proc.time()
      
      # Record time elapsed
      time <- round((t2 - t1)[3]/60, 2)
      fwrite(
        data.frame(time = time, method = "scmet", n_cores = n_cores),  
        paste0("data/interim/sim_studies/benchmark_real_chr/scmet/output/modelTime_", 
               bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
               NV, "VMRs_sparseLevel", sparseLevel, 
               "_seed", seed, ".txt")
      )
      
      # # plot mean-overdispersion
      # gg1 <- scmet_plot_mean_var(obj = fit_obj, y = "gamma", 
      #                            task = NULL, show_fit = TRUE)
      # gg2 <- scmet_plot_mean_var(obj = fit_obj, y = "epsilon", 
      #                            task = NULL, show_fit = TRUE)
      # cowplot::plot_grid(gg1, gg2, ncol = 2)
      # ggsave(paste0("plots/sim_studies/benchmark_real_chr/scmet/overdispersion_vs_mean_", 
      #               bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
      #               NV, "VMRs_sparseLevel", sparseLevel, 
      #               "_seed", seed, ".png"), width = 10, height = 5)
      
      # # store results from various efdr level
      # for (efdr in c(0.01, 0.02, 0.05, seq(0.1,0.9,0.1), 0.99)) {
      #   
      #   fit_obj <- scmet_hvf(scmet_obj = fit_obj, delta_e = 1-efdr, efdr = efdr)
      #   
      #   fwrite(
      #     fit_obj$hvf$summary, 
      #     paste0(
      #       "data/interim/sim_studies/benchmark_real_chr/scmet/output/pseudoChr_",
      #       subtype, "_", chromosome, "_", 
      #       bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
      #       NV, "VMRs_sparseLevel", sparseLevel, 
      #       "_seed", seed, "_", efdr, "efdr_scmetOutput.csv"
      #     ), 
      #     row.names = FALSE, col.names = TRUE
      #   )
      # }
    }
  }
}
