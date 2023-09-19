source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
library(argparse)
n_cores <- 22
register(MulticoreParam(workers = n_cores))

p <- ArgumentParser(description = 'Run vmrseq on modified real chromosome')
p$add_argument('--N',           type = "integer",      help = 'Number of cells')
p$add_argument('--NP',          type = "integer",      help = 'Number of subpopulation')
p$add_argument('--sparseLevel', type = "integer",      help = 'Sparse level')
p$add_argument('--alpha',       type = "double",       help = 'Alpha level')
args <- p$parse_args(commandArgs(TRUE))
N           <- args$N
NP          <- args$NP
sparseLevel <- args$sparseLevel
alpha       <- as.numeric(args$alpha)

cat("N =", N, "NP =", NP, "\n")
cat("alpha = ")

NV <- 2000
subtype <- "IT-L23_Cux1"
chromosome <- "chr1"
seed <- 2022

# load input
dir <- paste0("data/interim/sim_studies/benchmark_sim_chr/simulated/simChr_",
              subtype, "_", chromosome, "_", 
              N, "cells_", NP, "subpops_", 
              NV, "VMRs_sparseLevel", sparseLevel, 
              "_seed", seed)
SE <- loadHDF5SummarizedExperiment(dir)

# QC: remove sites with across-cell coverage < 3
total <- rowSums(assays(SE)[[1]]>=0, na.rm = T)
SE <- subset(SE, total >= 3)

# run model
gr <- vmrseq.smooth(SE)
t1 <- proc.time()
fit <- vmrseq.fit(gr, alpha)
t2 <- proc.time()

# record time elapsed
time <- round((t2 - t1)[3]/60, 2)
fwrite(
  data.frame(time = time, method = "vmrseq", n_cores = n_cores),
  paste0("data/interim/sim_studies/benchmark_sim_chr/vmrseq/output/modelTime_",
         N, "cells_", NP, "subpops_",
         NV, "VMRs_sparseLevel", sparseLevel, 
         "_alpha", alpha, ".txt")
)

# save model output
saveRDS(fit, paste0(
  "data/interim/sim_studies/benchmark_sim_chr/vmrseq/output/simChr_",
  subtype, "_", chromosome, "_",
  N, "cells_", NP, "subpops_",
  NV, "VMRs_sparseLevel", sparseLevel, 
  "_alpha", alpha, "_seed", seed, "_vmrseqOutput.rds"
))

