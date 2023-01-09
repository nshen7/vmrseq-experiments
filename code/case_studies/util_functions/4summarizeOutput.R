.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(HDF5Array)
library(SummarizedExperiment)

# ==== main ====

### Concatenate GRanges of all CpG sites
summarizeAllSites <- function(read_dir, write_dir) {
  se_dirs <- paste0(read_dir, "vmrseq/input/chr", 1:19)
  sites.gr <- do.call(
    c, 
    lapply(se_dirs, function(dir) loadHDF5SummarizedExperiment(dir) %>% granges())
  )
  saveRDS(sites.gr, paste0(write_dir, "cpg_sites.rds"))
}

### Summarize single-site methylation
summarizeOutputSite <- function(read_dir, write_dir, methods) {
  
  se_dirs <- paste0(read_dir, "vmrseq/input/chr", 1:19)
  
  ### vmrseq
  if ("vmrseq" %in% methods) {
    # dirs_vseq <- paste0(read_dir, "vmrseq/output/chr", 1:19, ".rds")
    dirs_vseq <- paste0(read_dir, "vmrseq/output/chr", 1:19, ".rds")
    
    # Concatenate results from 19 chromosomes
    res_vseq <- list(
      vmr.ranges = do.call(c, lapply(dirs_vseq, function(dir) readRDS(dir)$vmr.ranges)),
      cr.ranges = do.call(c, lapply(dirs_vseq, function(dir) readRDS(dir)$cr.ranges))
    )
    
    # Subset sites in detected regions
    ss_vseq1.se <- getSiteSummary(
      se_dirs = se_dirs, 
      regions.gr = res_vseq$vmr.ranges
    )
    saveHDF5SummarizedExperiment(ss_vseq1.se, paste0(write_dir, "vmrseq_siteSummary_sparseRep_vmrs"), replace=TRUE)
    ss_vseq2.se <- getSiteSummary(
      se_dirs = se_dirs, 
      regions.gr = res_vseq$cr.ranges
    )
    saveHDF5SummarizedExperiment(ss_vseq2.se, paste0(write_dir, "vmrseq_siteSummary_sparseRep_crs"), replace=TRUE)
    print("Finihsed vmrseq.")
  }
  
  ### scbs
  if ("scbs" %in% methods) {
    # Convert scbs output to GRanges obj
    res_scbs <- fread(paste0(read_dir, "scbs/output/scbs_output_0.02vt.bed"))
    res_scbs.gr <- with(
      res_scbs,
      GRanges(seqnames = V1, ranges = IRanges(start = V2, end = V3), mcols = data.frame("var" = V4))
    )
    
    ss_scbs.se <- getSiteSummary(
      se_dirs = se_dirs, 
      regions.gr = res_scbs.gr
    )
    saveHDF5SummarizedExperiment(ss_scbs.se, paste0(write_dir, "scbs_siteSummary_sparseRep_vmrs"), replace=TRUE)
    print("Finihsed scbs.")
  }
  
  ### smallwood 
  if ("smallwood" %in% methods) {
    all_smwd.gr <- readRDS(paste0(read_dir, "smallwood/output/smallwood_output_varLowerBound.rds"))
    
    # Take top 2% windows with highest lower bound and merge them as VMRs
    n <- round(0.02 * length(all_smwd.gr))
    ind <- order(all_smwd.gr$var_lb, decreasing = TRUE)[1:n]
    res_smwd.gr <- all_smwd.gr[ind]
    
    ss_smwd.se <- getSiteSummary(
      se_dirs = se_dirs,
      regions.gr = res_smwd.gr %>% GenomicRanges::reduce()
    )
    saveHDF5SummarizedExperiment(ss_smwd.se, paste0(write_dir, "smallwood_siteSummary_sparseRep_vmrs"), replace=TRUE)
    print("Finihsed smallwood.")
  }
  
  ### scMET
  if ("scmet" %in% methods) {
    # Read in HVF detecetd by scmet
    hvf_scmet <- fread(paste0(read_dir, "scmet/output/scmet_summary_output_efdr0.1.csv")) %>% 
      filter(is_variable) %>%
      mutate(chr = gsub("chr(.*)_window_.*", replacement = "\\1", x = feature_name) %>% as.integer()) %>%
      mutate(feat_index = gsub("chr.*_window_(.*)", replacement = "\\1", x = feature_name) %>% as.integer())
    
    # Convert into GRanges object
    feat_md_scmet <- readRDS(paste0(read_dir, "scmet/input/scmet_feature_metadata.rds"))
    hvf_scmet.gr <- do.call(c, lapply(
      1:nrow(hvf_scmet), 
      function(i) feat_md_scmet[[hvf_scmet$chr[i]]][hvf_scmet$feat_index[i]]
    ))
    
    # Add metadata info to GRanges object
    values(hvf_scmet.gr) <- hvf_scmet[,c("mu", "gamma", "epsilon", "tail_prob")]
    seqlevels(hvf_scmet.gr) <- paste0("chr", 1:19)
    hvf_scmet.gr <- sort(hvf_scmet.gr)
    
    ss_scmet.se <- getSiteSummary(
      se_dirs = se_dirs,
      regions.gr = hvf_scmet.gr
    )
    saveHDF5SummarizedExperiment(ss_scmet.se, paste0(write_dir, "scmet_siteSummary_sparseRep_vmrs"), replace=TRUE)
  }
}


### Summarize region methylation
summarizeOutputRegion <- function(read_dir, write_dir, methods) {
  
  se_dirs <- paste0(read_dir, "vmrseq/input/chr", 1:19)
  
  ### vmrseq
  if ("vmrseq" %in% methods) {
    dirs_vseq <- paste0(read_dir, "vmrseq/output/chr", 1:19, ".rds")
    
    # Concatenate results from 19 chromosomes
    res_vseq <- list(
      vmr.ranges = do.call(c, lapply(dirs_vseq, function(dir) readRDS(dir)$vmr.ranges)),
      cr.ranges = do.call(c, lapply(dirs_vseq, function(dir) readRDS(dir)$cr.ranges))
    )

    rs_vseq1.se <- getRegionSummary(
      se_dirs = se_dirs,
      regions.gr = res_vseq$vmr.ranges
    )
    saveHDF5SummarizedExperiment(rs_vseq1.se, paste0(write_dir, "vmrseq_regionSummary_vmrs"), replace=TRUE)
    rs_vseq2.se <- getRegionSummary(
      se_dirs = se_dirs, 
      regions.gr = res_vseq$cr.ranges
    )
    saveHDF5SummarizedExperiment(rs_vseq2.se, paste0(write_dir, "vmrseq_regionSummary_crs"), replace=TRUE)
  }
  
  ### scbs
  if ("scbs" %in% methods) {
    # Convert scbs output to GRanges obj
    res_scbs <- fread(paste0(read_dir, "scbs/output/scbs_output_0.02vt.bed"))
    res_scbs.gr <- with(
      res_scbs,
      GRanges(seqnames = V1, ranges = IRanges(start = V2, end = V3), mcols = data.frame("var" = V4))
    )
    
    rs_scbs.se <- getRegionSummary(
      se_dirs = se_dirs, 
      regions.gr = res_scbs.gr
    )
    saveHDF5SummarizedExperiment(rs_scbs.se, paste0(write_dir, "scbs_regionSummary_vmrs"), replace=TRUE)
  }
  
  ### smallwood 
  if ("smallwood" %in% methods) {
    all_smwd.gr <- readRDS(paste0(read_dir, "smallwood/output/smallwood_output_varLowerBound.rds"))
    
    # Take top 2% windows with highest lower bound (they might overlap)
    n <- round(0.02 * length(all_smwd.gr))
    ind <- order(all_smwd.gr$var_lb, decreasing = TRUE)[1:n]
    res_smwd.gr <- all_smwd.gr[ind] 
    
    rs_smwd.se <- getRegionSummary(
      se_dirs = se_dirs,
      regions.gr = res_smwd.gr
    )
    saveHDF5SummarizedExperiment(rs_smwd.se, paste0(write_dir, "smallwood_regionSummary_vmrs"), replace=TRUE)
  }
  
  ### scMET
  if ("scmet" %in% methods) {
    # Read in HVF detecetd by scmet
    hvf_scmet <- fread(paste0(read_dir, "scmet/output/scmet_summary_output_efdr0.1.csv")) %>% 
      filter(is_variable) %>%
      mutate(chr = gsub("chr(.*)_window_.*", replacement = "\\1", x = feature_name) %>% as.integer()) %>%
      mutate(feat_index = gsub("chr.*_window_(.*)", replacement = "\\1", x = feature_name) %>% as.integer())

    # Convert into GRanges object
    feat_md_scmet <- readRDS(paste0(read_dir, "scmet/input/scmet_feature_metadata.rds"))
    hvf_scmet.gr <- do.call(c, lapply(
      1:nrow(hvf_scmet), 
      function(i) feat_md_scmet[[hvf_scmet$chr[i]]][hvf_scmet$feat_index[i]]
    ))
     
    # Add metadata info to GRanges object
    values(hvf_scmet.gr) <- hvf_scmet[,c("mu", "gamma", "epsilon", "tail_prob")]
    seqlevels(hvf_scmet.gr) <- paste0("chr", 1:19)
    hvf_scmet.gr <- sort(hvf_scmet.gr)
    
    rs_scmet.se <- getRegionSummary(
      se_dirs = se_dirs,
      regions.gr = hvf_scmet.gr
    )
    saveHDF5SummarizedExperiment(rs_scmet.se, paste0(write_dir, "scmet_regionSummary_vmrs"), replace=TRUE)
  }
}


preparePDclustInput <- function(se, sparse) {
  if (sparse) {
    
  } else {
    break
  }
}




# ==== utils ====
getSiteSummary <- function(se_dirs, regions.gr) {
  
  if (any(duplicated(se_dirs))) stop("Duplicates exists in `se_dirs`.")
  
  for (se_dir in se_dirs) {
    se <- loadHDF5SummarizedExperiment(dir = se_dir)
    hits <- findOverlaps(granges(se), regions.gr)
    idx <- queryHits(hits)
    
    temp_sites.se <- se[idx, ]
    temp_sites.gr <- granges(temp_sites.se)
    temp_M_mat <- assays(temp_sites.se)$M_mat %>% as("sparseMatrix")
    if (se_dir == se_dirs[1]) {
      sites.gr <- temp_sites.gr
      M_mat <- temp_M_mat
    } else {
      sites.gr <- c(sites.gr, temp_sites.gr)
      M_mat <- rbind(M_mat, temp_M_mat)
    }
    print(se_dir)
    print(length(temp_sites.se))
  }
  
  sites.se <- SummarizedExperiment(assays = list("M_mat" = M_mat))
  return(sites.se)
}


getRegionSummary <- function(se_dirs, regions.gr) {
  
  if (any(duplicated(se_dirs))) stop("Duplicates exists in `se_dirs`.")
  
  for (se_dir in se_dirs) {
    
    se <- loadHDF5SummarizedExperiment(dir = se_dir)
    hits <- findOverlaps(granges(se), regions.gr)
    idx <- unique(subjectHits(hits))
    
    computeRegionStats <- function(i, type) { # i th feature/window
      mat <- assays(se)$M_mat[queryHits(hits)[subjectHits(hits)==i], ] %>% as("sparseMatrix")
      if (type == "M") return(round(colSums(mat))) 
      else if (type == "Cov") return(colSums(mat > 0))
      else stop("Wrong 'type' value. Either 'Cov' or 'M'.")
    }
    M <- do.call(
      rbind,
      bplapply(idx, computeRegionStats, type = "M")
    )
    Cov <- do.call(
      rbind,
      bplapply(idx, computeRegionStats, type = "Cov")
    )
    
    temp_regions.se <- SummarizedExperiment(
      assays = list("M" = M, "Cov" = Cov),
      rowRanges = regions.gr[idx]
    )
    
    if (se_dir == se_dirs[1]) regions.se <- temp_regions.se else regions.se <- rbind(regions.se, temp_regions.se)
    cat("Finished", se_dir, "\n")
  }
  return(regions.se)
}
