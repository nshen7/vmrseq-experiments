library(data.table)
library(readxl)
library(tidyverse)
setwd(here::here())
# source("code/liu2021_data_processing/helper_functions.R")

### import and clean metadata
md0 <- as.data.table(read_xlsx("data/metadata/metadata_liu2021/Liu2021_cell_metadata.xlsx", skip = 15))
colnames(md0)[1] <- "CellID"

sample_geo <- fread("data/metadata/metadata_liu2021/GSE132489_GEOaccession_with_CEMBA_ID.txt") %>%
  dplyr::rename(CEMBA_ID = `CEMBA ID`, GEO_accession = `GEO accession`) %>%
  mutate(Sample = gsub(pattern = "CEMBA(.*)_(.*)", replacement = "\\2_\\1", x = CEMBA_ID)) 
md <- merge.data.table(md0, sample_geo, by = "Sample", all = T)

# > colnames(md)
# [1] "Sample"          "CellID"          "CCC_Frac"        "CG_Frac"        
# [5] "CH_Frac"         "MappingRate"     "FinalReads"      "GeneCov>20"     
# [9] "Chrom100kCov>20" "Pass QC"         "RegionName"      "MajorRegion"    
# [13] "SubRegion"       "CellClass"       "MajorType"       "SubType"        
# [17] "L1UMAP_0"        "L1UMAP_1"        "L1TSNE_0"        "L1TSNE_1"       
# [21] "L2UMAP_0"        "L2UMAP_1"        "L2TSNE_0"        "L2TSNE_1"       
# [25] "L3UMAP_0"        "L3UMAP_1"        "L3TSNE_0"        "L3TSNE_1"       
# [29] "GEO_accession"   "CEMBA_ID"  

cellid_to_file <- as.data.table(read.csv("data/metadata/metadata_liu2021/cell_id_to_allc_name.csv", header = F))
colnames(cellid_to_file) <- c("CellID", "fullCEMBA")
cellid_to_file$subCEMBA <- gsub(pattern = "-", replacement = "_", x = cellid_to_file$fullCEMBA)
cellid_to_file$subCEMBA <- gsub(pattern = "ad", replacement = "AD", x = cellid_to_file$subCEMBA)
cellid_to_file$subCEMBA <- gsub(pattern = ".*(CEMBA[^_].*CEMBA[^_].*AD\\d{3}).*", replacement = "\\1", x = cellid_to_file$subCEMBA)
metadata <- merge.data.table(cellid_to_file, md, by = "CellID")

write.csv(metadata, "data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv",
          row.names = FALSE, quote = FALSE)
