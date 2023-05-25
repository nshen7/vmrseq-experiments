source("code/SETPATHS.R")
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(BiocParallel)
# register(MulticoreParam(workers = n_cores))

formatDataScmet <- function(read_dir, write_dir, chrs, bp_size, n_cores) {
  
  if (!file.exists(write_dir)) dir.create(write_dir)
  
  all_wds.gr <- NULL
  for (chr in chrs) {
    
    se <- loadHDF5SummarizedExperiment(dir = paste0(read_dir, chr))
    message(1)
    
    # Format into scMET input
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
      feat.df <- data.frame(
        Feature = paste0(chr, "_window_", i),
        Cell = paste0("cell_",1:length(total_reads)),
        total_reads = total_reads,
        met_reads = met_reads
      ) %>% dplyr::filter(total_reads > 0)
      return(feat.df)
    }
    
    feats.df <- do.call(
      rbind,
      parallel::mclapply(unique(subjectHits(hits)), computeFeature, mc.cores = n_cores)
    )
    message(3)
    
    # Save input for scMET 
    append <- ifelse(chr != chrs[1], yes = TRUE, no = FALSE) 
    # print(head(feats.df)); print(str(feats.df)); print(dim(feats.df))
    fwrite(feats.df, file = paste0(write_dir, "/scmet_input_Y.txt.gz"), 
           col.names = F, quote = F, append = append)
    print(4)
    
    all_wds.gr <- c(all_wds.gr, wds.gr)
    message(5)
    
    cat(chr, " ")
  }
  
  # Save feature metadata
  saveRDS(all_wds.gr, file = paste0(write_dir, "/scmet_feature_metadata.rds"))
}
