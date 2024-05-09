suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(parallel))

source("code/SETPATHS.R")
source("code/raw_data_process/liu2021_data_processing/data_processing_subtypeAcrossSample_helper_functions.R")

## import metadata
metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv") %>%
  filter(!is.na(GEO_accession) & !is.na(SubType) & !is.na(FilePath))
subtype_smr <- metadata[, .(.N), by = .(CellClass, SubType)] %>% 
  arrange(desc(N)) %>% 
  filter(!grepl("Outlier", SubType)) %>%
  mutate(Bin = cut(N, breaks = seq(0, ceiling(max(N)/100)*100, 100), labels = F)) %>%
  group_by(Bin) %>% mutate(nSubTypeInBin = n())

hist(subtype_smr$N, breaks = seq(0, ceiling(max(subtype_smr$N)/100)*100, 100))

## Select (at most) 3 subtypes from each bin of 100 in N cell
set.seed(2022)
subtype_smr_sel <- rbind(subtype_smr %>% 
                           filter(nSubTypeInBin <= 3),
                         subtype_smr %>% 
                           filter(nSubTypeInBin > 3) %>%
                           group_by(Bin) %>% 
                           sample_n(3) %>% 
                           mutate(nSubTypeInBin = n())
                         ) %>% 
  arrange(desc(N))

# nrow(subtype_smr_sel) # = 53 subtypes

# > subtype_smr_sel$SubType
# [1] "IT-L23 Cux1"      "DG dg-all"        "CT-L6 Il1rap"     "OLF-Exc Bmpr1b"  
# [5] "CT-L6 Megf9"      "CA1 Chrm3"        "IT-L4 Astn2"      "IT-L6 Man1c1"    
# [9] "IT-L5 Cdh8"       "OLF-Exc Sgcd"     "IT-L23 Tenm2"     "CA3 Cadm2"       
# [13] "IT-L5 Etv1"       "ODC odc-large"    "MSN-D2 Slc24a2"   "PAL-Inh Meis2"   
# [17] "IT-L5 Grik3"      "OLF-Exc Lrrtm3"   "IT-L4 Shc3"       "IT-L23 Foxp1"    
# [21] "IT-L23 Ptprt"     "MSN-D1 Hrh1"      "IT-L6 Fstl4"      "ASC str-hpf"     
# [25] "OLF Trpc4"        "OLF-Exc Unc13c"   "IT-L6 Cadps2"     "ODC odc-small"   
# [29] "MGC mgc-all"      "CGE-Lamp5 Grk5"   "OLF-Exc Cdh9"     "ASC mid"         
# [33] "MGE-Sst Frmd6"    "OPC opc-large"    "MGE-Pvalb Thsd7a" "MSN-D2 Col14a1"  
# [37] "MSN-D1 Khdrbs3"   "CLA Cdh8"         "L6b Nrp2"         "Foxp2 Inpp4b"    
# [41] "ASC cortex-olf"   "PT-L5 Ptprt"      "MGE-Sst Dock4"    "OLF-Exc Rmst"    
# [45] "MGE-Sst Bmper"    "D1L-PAL Plcxd3"   "ANP anp-olf-cnu"  "PAL-Inh Onecut2" 
# [49] "PC pc-all"        "MGE-Sst Unc5b"    "OPC opc-small"    "OLF Mapk10"      
# [53] "OLF Xkr6"  


# ==== run ====

for (i in 5) wrapper(subtype_smr_sel$SubType[i], mc.cores = 10) # 3826303.pbsha.ib.sockeye
