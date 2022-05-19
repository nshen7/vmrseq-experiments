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


# cell_counts_inh <- cell_counts %>%
#   filter(CellClass=="Inh") %>%
#   arrange(desc(N)) %>%
#   filter(!is.na(SubType), N >= 100)
# cell_counts_inh[which(!duplicated(SubType)&!duplicated(Sample))]
#       Sample CellClass           SubType GEO_accession   N 
# 1: 5H_181015       Inh     PAL-Inh Meis2     GSE134806 498 (processed) (qced)
# 2: 1C_180208       Inh         OLF Trpc4     GSE131836 439 (processed) (qced)
# 3: 5E_180925       Inh       MSN-D1 Hrh1     GSE134653 369 (processed) (qced)
# 4: 4D_171219       Inh    MSN-D2 Slc24a2     GSE131406 353 (processed) (qced)
# 5: 5J_190207       Inh    LSX-Inh Dock10     GSE135307 241 (processed) (qced)
# 6: 4H_180911       Inh     PAL-Inh Ptprd     GSE134617 164 (processed) (qced)
# 7: 5F_181220       Inh     MSN-D1 Plxnc1     GSE135169 158 (processed) (qced) (not all cells are available in raw data)
# 8: 4F_180329       Inh D1L-Fstl4 Sipa1l2     GSE132445 117 (processed) (qced)

# read_dirs <- c("data/processed/processed_liu2021/sample5J_190207_GSE135307_subtype_LSX-Inh_Dock10_241cells",
#                "data/processed/processed_liu2021/sample4H_180911_GSE134617_subtype_PAL-Inh_Ptprd_164cells",
#                "data/processed/processed_liu2021/sample5F_181220_GSE135169_subtype_MSN-D1_Plxnc1_123cells",
#                "data/processed/processed_liu2021/sample4F_180329_GSE132445_subtype_D1L-Fstl4_Sipa1l2_117cells")

for (read_dir in read_dirs) {
  write_dir <- paste0(read_dir, "_qced")
  saveHDF5SummarizedExperiment(qcOnSE(read_dir), dir = write_dir, replace = T)
}
