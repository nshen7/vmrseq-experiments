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
NV <- 2000

for (N in c(2000)) {
  # for (NP in c(2)) {
  # for (NP in c(3)) {
  # for (NP in c(4)) {
  # for (NP in c(5)) {
  for (NP in c(8)) {
  # for (NP in c(12)) {
  # for (NP in c(20)) {
    for (sparseLevel in 3) {
    # for (sparseLevel in 1:3) {
      
      # Load raw data
      load_dir <- paste0(
        "data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
        subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
        NV, "VMRs_sparseLevel", sparseLevel, 
        "_seed", seed
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
        M_mat <- assays(se)$M_mat[queryHits(hits)[subjectHits(hits)==i],]
        met_reads <- colSums(M_mat)
        total_reads <- colSums(M_mat > 0)
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
        NV, "VMRs_sparseLevel", sparseLevel, 
        "_seed", seed, ".rds"
      )
      saveRDS(wds.gr, file = wds_dir)
      # Save input for scMET
      write_dir <- paste0(
        input_folder,
        "/pseudoChr_", subtype, "_", chromosome, "_",
        bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
        NV, "VMRs_sparseLevel", sparseLevel, 
        "_seed", seed, ".txt"
      )
      fwrite(feats.df, file = write_dir, quote = F)
      R.utils::gzip(write_dir, remove = T, overwrite = T)
      
      cat("N =", N, "; NP =", NP, "\n")
    }
  }
}
