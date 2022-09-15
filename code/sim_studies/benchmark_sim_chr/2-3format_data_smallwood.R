.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") 
library(tidyverse)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)

register(MulticoreParam(workers = 30))

subtype <- "IT-L23_Cux1"
chromosome <- "chr1"
seed <- 2022
bp_window <- 3000
bp_step <- 600
NV <- 2000

# for (N in c(100, 500, 1000, 2000)) {
# for (N in c(100, 500, 1000)) {
for (N in c(200)) {
  # for (NP in c(2,3,4,5,8,12,20)) {
  for (NP in c(4)) {
    
    # Load raw data
    load_dir <- paste0(
      "data/interim/sim_studies/benchmark_sim_chr/simulated/simChr_",
      subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed
    )
    se <- loadHDF5SummarizedExperiment(dir = load_dir)
    
    # Format into scMET input
    cuts <- seq(start(se)[1], start(se)[length(se)], bp_step)
    len <- length(cuts)
    cuts <- cuts[-((len-3):len)]
    wds.gr <- GRanges(
      seqnames = seqnames(se)[1],
      ranges = IRanges(start = cuts,
                       end = c(cuts[-length(cuts)]+bp_window, start(se)[length(se)]))
    )
    hits <- findOverlaps(granges(se), wds.gr)

    getFeature <- function(i, type) { # i th feature/window
      feat.se <- se[queryHits(hits)[subjectHits(hits)==i],]
      if (type == "Cov") {
        return(colSums(assays(feat.se)$M_mat >= 0, na.rm = T)) 
      } else if (type == "M") {
        return(colSums(assays(feat.se)$M_mat, na.rm = T))
      } else {
        stop("Wrong 'type' value. Either 'Cov' or 'M'.")
      }
    }
    
    M <- do.call(
      rbind,
      bplapply(unique(subjectHits(hits)), getFeature, type = "M")
      # bplapply(1:10000, getFeature, type = "M")
    )
    Cov <- do.call(
      rbind,
      bplapply(unique(subjectHits(hits)), getFeature, type = "Cov")
      # bplapply(1:10000, getFeature, type = "Cov")
    )
    
    
    input_folder <- paste0("data/interim/sim_studies/benchmark_sim_chr/smallwood/input")
    write_dir <- paste0(
      input_folder,
      "/simChr_", subtype, "_", chromosome, "_",
      N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed
    )
    
    saveHDF5SummarizedExperiment(
      x = SummarizedExperiment(
        assays = list("M" = M, "Cov" = Cov),
        rowRanges = wds.gr[unique(subjectHits(hits))]
        # rowRanges = wds.gr[1:10000]
      ), 
      dir = write_dir, 
      replace = TRUE
    )
    cat("N =", N, "; NP =", NP, "\n")
  }
}
