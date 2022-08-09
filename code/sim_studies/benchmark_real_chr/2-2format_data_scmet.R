.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") 
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(BiocParallel)

# ==== util ====
formatForScmet <- function(se, bp_size) {
  
  cuts <- seq(start(se)[1], start(se)[length(se)], bp_size)
  wds.gr <- GRanges(
    seqnames = seqnames(se)[1],
    ranges = IRanges(start = cuts,
                     end = c(cuts[-1]-1, start(se)[length(se)]))
  )
  hits <- findOverlaps(granges(se), wds.gr)
  
  computeFeature <- function(i) { # i th feature/window
    feat.se <- se[queryHits(hits)[subjectHits(hits)==i],]
    total_reads <- colSums(assays(feat.se)$MF >= 0, na.rm = T)
    met_reads <- colSums(assays(feat.se)$MF, na.rm = T)
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
    bplapply(1:length(wds.gr), computeFeature)
  )
  return(feats.df)
}

# ==== main ====
bp_size <- 20000
for (N in c(100, 500, 1000, 5000)) {
  for (NP in c(2,3,4,5,8,12,20)) {

    # Load raw data
    load_dir <- paste0(
      "data/interim/sim_studies/benchmark_pseudo_chr/raw/pseudo_chr_",
      N, "cells_",
      NP, "subpops"
    )
    se <- loadHDF5SummarizedExperiment(dir = load_dir)
    
    # Format into scMET input
    feats.df <- formatForScmet(se, bp_size = bp_size)
    
    # Write file
    input_folder <- paste0("data/interim/sim_studies/benchmark_pseudo_chr/scmet/input")
    write_dir <- paste0(input_folder, "/scmet_input_", bp_size/1000, "kbWindow_",
                        N, "cells_", NP, "subpops.txt")
    fwrite(feats.df, file = write_dir, quote = F)
    R.utils::gzip(write_dir, remove = T, overwrite = T)
    
    cat("N =", N, "; NP =", NP, "\n")
  }
}
