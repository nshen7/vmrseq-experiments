.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") 
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(BiocParallel)


subtype <- "IT-L23_Cux1"
chromosome <- "chr1"
seed <- 2022
bp_size <- 20000
NV <- 1500

# for (N in c(100, 500, 1000, 2000)) {
# for (N in c(100, 500, 1000)) {
for (N in c(200)) {
  # for (NP in c(2,3,4,5,8,12,20)) {
  for (NP in c(4)) {
    
    # Load raw data
    load_dir <- paste0(
      "data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
      subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed
    )
    se <- loadHDF5SummarizedExperiment(dir = load_dir)
    
    # Format into scMET input
    cuts <- seq(start(se)[1], start(se)[length(se)], bp_size)
    wds.gr <- GRanges(
      seqnames = seqnames(se)[1],
      ranges = IRanges(start = cuts,
                       end = c(cuts[-1]-1, start(se)[length(se)]))
    )
    hits <- findOverlaps(granges(se), wds.gr)
    
    computeFeature <- function(i) { # i th feature/window
      feat.se <- se[queryHits(hits)[subjectHits(hits)==i],]
      total_reads <- colSums(assays(feat.se)$M_mat >= 0, na.rm = T)
      met_reads <- colSums(assays(feat.se)$M_mat, na.rm = T)
      feat.df <- data.frame(
        Feature = paste0("Window_", i),
        Cell = paste0("Cell_",1:length(total_reads)),
        total_reads = total_reads,
        met_reads = met_reads
      ) %>% dplyr::filter(total_reads > 0)
      return(feat.df)
    }
    
    feats.df <- do.call(
      rbind,
      bplapply(unique(subjectHits(hits)), computeFeature)
    )
    
    input_folder <- paste0("data/interim/sim_studies/benchmark_real_chr/scmet/input")
    # Save feature metadata
    wds_dir <- paste0(
      input_folder, 
      "/features_", subtype, "_", chromosome, "_", 
      bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed, ".rds"
    )
    saveRDS(wds.gr, file = wds_dir)
    # Save input for scMET
    write_dir <- paste0(
      input_folder,
      "/pseudoChr_", subtype, "_", chromosome, "_",
      bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed, ".txt"
    )
    fwrite(feats.df, file = write_dir, quote = F)
    R.utils::gzip(write_dir, remove = T, overwrite = T)
    
    cat("N =", N, "; NP =", NP, "\n")
  }
}
