suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(DelayedArray))
suppressPackageStartupMessages(library(DelayedMatrixStats))
suppressPackageStartupMessages(library(tidyverse))
setAutoBlockSize(size=1e9)

setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")

qcOnSE <- function(se_dir, min_cell_cov = 5, rm_allel_var = T){
  cells0.se <- loadHDF5SummarizedExperiment(dir = se_dir)
  print("Finished reading in processed pre-QC SE object.")

  cell_cov <- rowSums(assays(cells0.se)$Cov > 0, na.rm = T)
  cells.se <- cells0.se[cell_cov >= min_cell_cov, ]
  print("Finished removing low-cell-coverage sites.")
  
  if (rm_allel_var) {
    MF <- assays(cells.se)$M / assays(cells.se)$Cov
    itmd <- rowAnys(MF > 0 & MF < 1, na.rm = T)
    cells.se <- cells.se[!itmd, ]
    print("Finished removing sites with possible allelic variation.")
  }

  values(cells.se)$cell_cov <- rowSums(assays(cells.se)$Cov > 0, na.rm = T)
  values(cells.se)$cell_meth <- rowSums(assays(cells.se)$M > 0, na.rm = T)
  values(cells.se)$cell_MF <- values(cells.se)$cell_meth/values(cells.se)$cell_cov
  print("Finished computing important statistics.")
  
  return(cells.se)
}



#### Perform QC ====
read_dir <- "data/processed/processed_luo2017_human/subtype_hDL-1_144cells" # can't have '/' at the end
write_dir <- paste0(read_dir, "_qced")
saveHDF5SummarizedExperiment(qcOnSE(read_dir), dir = write_dir, replace = T)

read_dir <- "data/processed/processed_luo2017_human/subtype_hL5-4_162cells"
write_dir <- paste0(read_dir, "_qced")
saveHDF5SummarizedExperiment(qcOnSE(read_dir), dir = write_dir, replace = T)

read_dir <- "data/processed/processed_luo2017_human/subtype_hNdnf_173cells"
write_dir <- paste0(read_dir, "_qced")
saveHDF5SummarizedExperiment(qcOnSE(read_dir), dir = write_dir, replace = T)

read_dir <- "data/processed/processed_luo2017_human/subtype_hPv-1_175cells"
write_dir <- paste0(read_dir, "_qced")
saveHDF5SummarizedExperiment(qcOnSE(read_dir), dir = write_dir, replace = T)

read_dir <- "data/processed/processed_luo2017_human/subtype_hSst-2_163cells"
write_dir <- paste0(read_dir, "_qced")
saveHDF5SummarizedExperiment(qcOnSE(read_dir), dir = write_dir, replace = T)

read_dir <- "data/processed/processed_luo2017_human/subtype_hL23_873cells"
write_dir <- paste0(read_dir, "_qced")
saveHDF5SummarizedExperiment(qcOnSE(read_dir), dir = write_dir, replace = T)
