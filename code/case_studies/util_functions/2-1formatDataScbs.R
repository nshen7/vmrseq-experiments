source("code/SETPATHS.R")
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(recommenderlab)


formatDataScbs <- function(read_dir, write_dir, chrs, n_cores) {
  for (chr in chrs) {
    se <- loadHDF5SummarizedExperiment(dir = paste0(read_dir, chr))
    
    formatChr <- function(i) { # ith cell
      chr.se <- se[, i]
      M_vec <- dropNA2matrix(as(assays(chr.se)$M_mat, "dgCMatrix")) %>% as.vector()
      chr.df <- data.frame(
        seqnames = seqnames(chr.se),
        start = start(chr.se),
        end = end(chr.se),
        MF = M_vec*100,
        meth = M_vec,
        unmeth = 1-M_vec
      ) %>% filter(!is.na(MF))

      append <- ifelse(chr != chrs[1], yes = TRUE, no = FALSE) 
      fwrite(chr.df, append = append,
             file = paste0(write_dir, "/cell_", i, ".cov"),
             row.names = F, col.names = F, sep = "\t", quote = F)
    }
    
    formatChr(1)
    parallel::mclapply(1:ncol(se), formatChr, mc.cores = n_cores)
  }
}