suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(tidyverse))

setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
source("code/liu2021_data_processing/data_processing_helper_functions.R")

### import metadata
# metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv")
# cell_counts <- metadata[, .(.N), by = .(Sample, CellClass, SubType, GEO_accession)]


# ==== Inhibitory cell class ====

# cell_counts_inh <- cell_counts %>%
#   filter(CellClass=="Inh") %>%
#   arrange(desc(N)) %>%
#   filter(!is.na(SubType), N >= 100)
# cell_counts_inh[which(!duplicated(SubType)&!duplicated(Sample))]
#       Sample CellClass           SubType GEO_accession   N 
# 1: 5H_181015       Inh     PAL-Inh Meis2     GSE134806 498 (processed)
# 2: 1C_180208       Inh         OLF Trpc4     GSE131836 439 (processed)
# 3: 5E_180925       Inh       MSN-D1 Hrh1     GSE134653 369 (processed)
# 4: 4D_171219       Inh    MSN-D2 Slc24a2     GSE131406 353 (processed)
# 5: 5J_190207       Inh    LSX-Inh Dock10     GSE135307 241 (processed)
# 6: 4H_180911       Inh     PAL-Inh Ptprd     GSE134617 164 (processed)
# 7: 5F_181220       Inh     MSN-D1 Plxnc1     GSE135169 158 (processed) (not all cells are available in raw data)
# 8: 4F_180329       Inh D1L-Fstl4 Sipa1l2     GSE132445 117 (processed)

# wrapper(sample_name = "5J_190207", geo_accession = "GSE135307", subtype = "LSX-Inh Dock10", mc.cores = 16)
# wrapper(sample_name = "4H_180911", geo_accession = "GSE134617", subtype = "PAL-Inh Ptprd", mc.cores = 16)

# wrapper(sample_name = "5F_181220", geo_accession = "GSE135169", subtype = "MSN-D1 Plxnc1", mc.cores = 16)
# wrapper(sample_name = "4F_180329", geo_accession = "GSE132445", subtype = "D1L-Fstl4 Sipa1l2", mc.cores = 16)
