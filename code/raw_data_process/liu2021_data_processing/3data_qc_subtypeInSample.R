#### !!remember that the input folder name cannot have `/` at the end!!
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(DelayedArray))
suppressPackageStartupMessages(library(DelayedMatrixStats))
suppressPackageStartupMessages(library(tidyverse))
setAutoBlockSize(size=1e9)

source("code/SETPATHS.R")

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

# ==== Ecitatory cell class ====
### import metadata
# metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv")
# cell_counts <- metadata[, .(.N), by = .(Sample, CellClass, SubType, GEO_accession)]

# cell_counts_exc <- cell_counts %>%
#   filter(CellClass=="Exc") %>%
#   arrange(desc(N)) %>%
#   filter(!is.na(SubType), N >= 100) 
# cell_counts_exc[c(5, which(!duplicated(SubType)&!duplicated(Sample)))]
#         Sample CellClass        SubType GEO_accession    N
# 1:  11F_190214       Exc      DG dg-all     GSE135543  871 (processed) (qced)
# 2:   9J_190212       Exc      DG dg-all     GSE135484 1154 (deprecated; hard to process and might induce influential point in prior fitting)
# 3:   9H_190219       Exc      CA1 Chrm3     GSE135570  685 (processed) (qced)
# 4:   2E_180220       Exc OLF-Exc Bmpr1b     GSE131904  564 (processed) (qced)
# 5:   1A_180226       Exc    IT-L23 Cux1     GSE130553  506 (processed) (qced)
# 6:   8E_190716       Exc      CA3 Cadm2     GSE158066  418 (processed) (qced)
# 7:   5G_181008       Exc   OLF-Exc Sgcd     GSE134773  410 (processed) (qced)
# 8:   1B_180213       Exc   IT-L23 Ptprt     GSE131889  369 (processed) (qced)
# 9:   6A_190108       Exc   IT-L23 Tenm2     GSE135179  267 (processed) (qced)
# 10: 11E_190214       Exc      CA1 Ptprg     GSE135534  237 (processed) (qced)
# 11:  5B_180529       Exc    IT-L4 Astn2     GSE132692  237 (processed) (qced)


#### Perform QC ====
# wrapper(sample_name = "1B_180213", geo_accession = "GSE131889", subtype = "IT-L23 Ptprt", mc.cores = 16)
# wrapper(sample_name = "6A_190108", geo_accession = "GSE135179", subtype = "IT-L23 Tenm2", mc.cores = 16)
# wrapper(sample_name = "11E_190214", geo_accession = "GSE135534", subtype = "CA1 Ptprg", mc.cores = 16)
# wrapper(sample_name = "5B_180529", geo_accession = "GSE132692", subtype = "IT-L4 Astn2", mc.cores = 16)
# 
# read_dir <- "data/processed/processed_liu2021/sample1B_180213_GSE131889_subtype_IT-L23_Ptprt_369cells"
# write_dir <- paste0(read_dir, "_qced")
# saveHDF5SummarizedExperiment(qcOnSE(read_dir), dir = write_dir, replace = T)
# 
# read_dir <- "data/processed/processed_liu2021/sample6A_190108_GSE135179_subtype_IT-L23_Tenm2_267cells"
# write_dir <- paste0(read_dir, "_qced")
# saveHDF5SummarizedExperiment(qcOnSE(read_dir), dir = write_dir, replace = T)
# 
# read_dir <- "data/processed/processed_liu2021/sample11E_190214_GSE135534_subtype_CA1_Ptprg_237cells"
# write_dir <- paste0(read_dir, "_qced")
# saveHDF5SummarizedExperiment(qcOnSE(read_dir), dir = write_dir, replace = T)
# 
# read_dir <- "data/processed/processed_liu2021/sample5B_180529_GSE132692_subtype_IT-L4_Astn2_237cells"
# write_dir <- paste0(read_dir, "_qced")
# saveHDF5SummarizedExperiment(qcOnSE(read_dir), dir = write_dir, replace = T)

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
