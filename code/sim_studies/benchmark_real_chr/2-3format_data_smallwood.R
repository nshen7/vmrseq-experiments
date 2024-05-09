source("code/SETPATHS.R")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)

register(MulticoreParam(workers = 22))

subtype <- "IT-L23_Cux1"
chromosome <- "chr1"
seed <- 2022
bp_window <- 3000
bp_step <- 600
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
      cuts <- seq(start(se)[1], start(se)[length(se)], bp_step)
      len <- length(cuts)
      cuts <- cuts[-((len-3):len)]
      wds.gr <- GRanges(
        seqnames = seqnames(se)[1],
        ranges = IRanges(
          start = cuts,
          end = c(cuts[-length(cuts)]+bp_window, start(se)[length(se)])
        )
      )
      hits <- findOverlaps(granges(se), wds.gr)
      Indexes <- lapply(
        unique(subjectHits(hits)), 
        function(i) queryHits(hits)[subjectHits(hits)==i]
      )
      print("Finished finding hits.")
      print(Sys.time())
      
      # Divide M_mat into groups to read into RAM and save time
      grp_size <- 50000
      n_grp <- length(Indexes) %/% grp_size
      if (length(Indexes) %% grp_size > 0) n_grp <- n_grp + 1
      
      M <- NULL; Cov <- NULL
      for (i in 1:n_grp) {
        grp_idx <- (grp_size*(i-1)+1):min(grp_size*i, length(Indexes))
        head <- Indexes[[grp_idx[1]]]
        tail <- Indexes[[grp_idx[length(grp_idx)]]]
        site_idx <- head[1]:tail[length(tail)]
        grp.se <- se[site_idx,]
        grp_mat <- assays(grp.se)$M_mat %>% as("sparseMatrix")
        
        getFeature <- function(j, type) { # i th feature/window
          mat <- matrix(grp_mat[Indexes[[j]]-site_idx[1]+1, ], ncol = N)
          if (type == "M") return(round(colSums(mat))) 
          else if (type == "Cov") return(colSums(mat > 0))
          else stop("Wrong 'type' value. Either 'Cov' or 'M'.")
        }
        M <- rbind(M, do.call(
          rbind,
          bplapply(grp_idx, getFeature, type = "M")
        ))
        Cov <- rbind(Cov, do.call(
          rbind,
          bplapply(grp_idx, getFeature, type = "Cov")
        ))
        cat(i, " ")
      }
      print("Finished computing features.")
      print(Sys.time())
      
      input_folder <- paste0("data/interim/sim_studies/benchmark_real_chr/smallwood/input")
      write_dir <- paste0(
        input_folder,
        "/pseudoChr_", subtype, "_", chromosome, "_",
        N, "cells_", NP, "subpops_",
        NV, "VMRs_sparseLevel", sparseLevel,
        "_seed", seed
      )

      saveHDF5SummarizedExperiment(
        x = SummarizedExperiment(
          assays = list("M" = M, "Cov" = Cov),
          rowRanges = wds.gr[unique(subjectHits(hits))]
        ),
        dir = write_dir,
        replace = TRUE
      )
      cat("N =", N, "; NP =", NP, "\n")
    }
  }
}
