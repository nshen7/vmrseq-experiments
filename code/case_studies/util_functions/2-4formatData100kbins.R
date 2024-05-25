source("code/SETPATHS.R")
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(BiocParallel)
# register(MulticoreParam(workers = n_cores))

formatData100kbins <- function(read_dir, write_dir, chrs, bp_size, n_cores) {
  
  if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)
  
  all_wds.gr <- NULL
  for (chr in chrs) {
    
    se <- loadHDF5SummarizedExperiment(dir = here(read_dir, chr))
    message(1)
    
    # Format into PCA input
    cuts <- seq(start(se)[1], start(se)[length(se)], bp_size)
    wds.gr <- GRanges(
      seqnames = seqnames(se)[1],
      ranges = IRanges(start = cuts,
                       end = c(cuts[-1]-1, start(se)[length(se)]))
    )
    message(2)
    hits <- findOverlaps(granges(se), wds.gr)
    
    computeFeature <- function(i) { # i th feature/window
      M_mat <- assays(se)$M_mat[queryHits(hits)[subjectHits(hits)==i], ]
      M_mat <- matrix(M_mat, ncol = ncol(se))
      met_reads <- as.integer(round(colSums(M_mat)))
      total_reads <- colSums(M_mat > 0)
      mf <- met_reads / total_reads
      return(mf)
    }
    
    mf.df <- do.call(
      rbind,
      parallel::mclapply(unique(subjectHits(hits)), computeFeature, mc.cores = n_cores)
    )
    message(3)
    
    # Save input 
    fwrite(mf.df, file = here(write_dir, "100kbins_input.txt.gz"), 
           col.names = F, quote = F, append = chr != chrs[1])
    print(4)
    
    all_wds.gr <- c(all_wds.gr, wds.gr)
    message(5)
    
    cat(chr, " ")
  }
  
  # Save feature metadata
  saveRDS(all_wds.gr, file = here(write_dir, "pca_feature_metadata.rds"))
}
