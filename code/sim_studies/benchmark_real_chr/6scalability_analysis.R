source("code/SETPATHS.R")

# ---- scalability analysis ----

NV <- 2000
alpha <- 0.05

time.df <- expand.grid(N = c(200, 500, 1000, 2000),
                       NP = c(2,3,4,5,8,12,20),
                       sparseLevel = 1:3,
                       time = 0)
for (i in 1:nrow(time.df)) {
  time <- fread(paste0("data/interim/sim_studies/benchmark_real_chr/vmrseq/output/modelTime_",
                       time.df$N[i], "cells_", time.df$NP[i], "subpops_",
                       NV, "VMRs_sparseLevel", time.df$sparseLevel[i],
                       "_alpha", alpha, ".txt"))
  
}

