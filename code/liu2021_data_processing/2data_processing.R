## !!! This script processes bed files and combine cells from same subtype and sample into SummarizedExperiment objects

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(tidyverse))

setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
source("code/liu2021_data_processing/data_processing_helper_functions.R")

### import metadata
# metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv")
# cell_counts <- metadata[, .(.N), by = .(Sample, CellClass, SubType, GEO_accession)]


# ==== Excitatory cell class ====
# cell_counts_exc <- cell_counts %>%
#   filter(CellClass=="Exc") %>%
#   arrange(desc(N)) %>%
#   filter(!is.na(SubType), N >= 100) 
# cell_counts_exc[c(5, which(!duplicated(SubType)&!duplicated(Sample)))]
#         Sample CellClass        SubType GEO_accession    N
# 1:  11F_190214       Exc      DG dg-all     GSE135543  871 (processed)
# 2:   9J_190212       Exc      DG dg-all     GSE135484 1154 (deprecated; hard to process and might induce influential point in prior fitting)
# 3:   9H_190219       Exc      CA1 Chrm3     GSE135570  685 (processed)
# 4:   2E_180220       Exc OLF-Exc Bmpr1b     GSE131904  564 (processed)
# 5:   1A_180226       Exc    IT-L23 Cux1     GSE130553  506 (processed)
# 6:   8E_190716       Exc      CA3 Cadm2     GSE158066  418 (processed)
# 7:   5G_181008       Exc   OLF-Exc Sgcd     GSE134773  410 (processed)
# 8:   1B_180213       Exc   IT-L23 Ptprt     GSE131889  369 (processed)
# 9:   6A_190108       Exc   IT-L23 Tenm2     GSE135179  267 (processed)
# 10: 11E_190214       Exc      CA1 Ptprg     GSE135534  237 (processed)
# 11:  5B_180529       Exc    IT-L4 Astn2     GSE132692  237 (processed)


# wrapper(sample_name = "1A_180226", geo_accession = "GSE130553", subtype = "IT-L23 Cux1", mc.cores = 16)
# wrapper(sample_name = "2E_180220", geo_accession = "GSE131904", subtype = "OLF-Exc Bmpr1b", mc.cores = 16)
# wrapper(sample_name = "8J_190716", geo_accession = "GSE158154", subtype = "DG dg-all", mc.cores = 16)
# wrapper(sample_name = "11F_190214", geo_accession = "GSE135543", subtype = "DG dg-all", mc.cores = 16)

# wrapper(sample_name = "8E_190716", geo_accession = "GSE158066", subtype = "CA3 Cadm2", mc.cores = 16)
# wrapper(sample_name = "5G_181008", geo_accession = "GSE134773", subtype = "OLF-Exc Sgcd", mc.cores = 16)
# wrapper(sample_name = "1B_180213", geo_accession = "GSE131889", subtype = "IT-L23 Ptprt", mc.cores = 16)
# wrapper(sample_name = "6A_190108", geo_accession = "GSE135179", subtype = "IT-L23 Tenm2", mc.cores = 16)
# wrapper(sample_name = "11E_190214", geo_accession = "GSE135534", subtype = "CA1 Ptprg", mc.cores = 16)
# wrapper(sample_name = "5B_180529", geo_accession = "GSE132692", subtype = "IT-L4 Astn2", mc.cores = 16)

### for examination of coverage vs. MF
# 4A_180205	Exc	IT-L23 Tenm2	GSE131766	190
# wrapper(sample_name = "4A_180205", geo_accession = "GSE131766", subtype = "IT-L23 Tenm2", mc.cores = 16)

