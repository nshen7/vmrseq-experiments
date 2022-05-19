## !!! This script processes raw bed files downloaded from GEO and extract out the CpGs and remove all non-CpGs.

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")

# ==== helper functions ====
subsetAndStrandFlip <- function(file_dir) {
  # cat("Start processing", cell_id, "\n")
  cell <- fread(
    file_dir,
    header = F,
    select = 1:6,
    colClasses = c("character", "integer", "character", "character",
                   "integer", "integer", "integer"))
  colnames(cell) <- c("chr","pos","strand","context","meth","total")
  
  ## we remove all non-CpG calls and only keeps autosomes.  
  cell <- cell[substr(context, start = 1, stop = 2) == "CG"  & chr%in%paste0("chr", 1:19)]
  cell$pos[cell$strand == "-"] <- cell$pos[cell$strand == "-"] - 1
  cell$strand <- "*"
  
  write_dir <- sub(".gz", "", file_dir)
  # write_dir <- sub("_RAW/","_RAW/Atest_", write_dir) ### for testing purpose
  write_tsv(cell, write_dir)
  R.utils::gzip(write_dir, remove = T, overwrite = T)
}


wrapper <- function(geo_accession, mc.cores){
  # cat("Start processing", cell_id, "\n")
  folder_dir <- paste0("data/raw_counts/raw_counts_Liu2021/",geo_accession,"_RAW/")
  file_names <- fread(paste0(folder_dir, "filelist.txt"))[-1,Name]
  file_dirs <- paste0(folder_dir, file_names)
  mclapply(file_dirs, subsetAndStrandFlip, mc.cores = mc.cores)
  # mclapply(file_dirs[1], subsetAndStrandFlip, mc.cores = mc.cores) ### for testing purpose
}



# ==== raw data subsetting ====

## import metadata
metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv")
cell_counts <- metadata[, .(.N), by = .(Sample, CellClass, SubType, GEO_accession)]

###### excitatory subtypes ######

cell_counts_exc <- cell_counts %>%
  filter(CellClass=="Exc") %>%
  arrange(desc(N)) %>%
  filter(!is.na(SubType), N >= 100)

(cell_counts_exc_selected <- cell_counts_exc[c(5, which(!duplicated(SubType)&!duplicated(Sample)))])
#         Sample CellClass        SubType GEO_accession    N
# 1:  11F_190214       Exc      DG dg-all     GSE135543  871 (raw downloaded) (CpG sites subsetted)
# 2:   9J_190212       Exc      DG dg-all     GSE135484 1154 (raw downloaded) (CpG sites subsetted)
# 3:   9H_190219       Exc      CA1 Chrm3     GSE135570  685 (raw downloaded) (CpG sites subsetted)
# 4:   2E_180220       Exc OLF-Exc Bmpr1b     GSE131904  564 (raw downloaded) (CpG sites subsetted)
# 5:   1A_180226       Exc    IT-L23 Cux1     GSE130553  506 (raw downloaded) (CpG sites subsetted)
# 6:   8E_190716       Exc      CA3 Cadm2     GSE158066  418 (raw downloaded) (CpG sites subsetted)
# 7:   5G_181008       Exc   OLF-Exc Sgcd     GSE134773  410 (raw downloaded) (CpG sites subsetted)
# 8:   1B_180213       Exc   IT-L23 Ptprt     GSE131889  369 (raw downloaded) (CpG sites subsetted)
# 9:   6A_190108       Exc   IT-L23 Tenm2     GSE135179  267 (raw downloaded) (CpG sites subsetted)
# 10: 11E_190214       Exc      CA1 Ptprg     GSE135534  237 (raw downloaded) (CpG sites subsetted)
# 11:  5B_180529       Exc    IT-L4 Astn2     GSE132692  237 (raw downloaded) (CpG sites subsetted)

# as.list(cell_counts_exc_selected$GEO_accession[6:11]) %>% map(wrapper, mc.cores = 16)


###### Inhibitory subtypes ######
cell_counts_inh <- cell_counts %>%
  filter(CellClass=="Inh") %>%
  arrange(desc(N)) %>%
  filter(!is.na(SubType), N >= 100)
(cell_counts_inh_selected <- cell_counts_inh[which(!duplicated(SubType)&!duplicated(Sample))])
#       Sample CellClass           SubType GEO_accession   N
# 1: 5H_181015       Inh     PAL-Inh Meis2     GSE134806 498 (raw downloaded) (CpG sites subsetted)
# 2: 1C_180208       Inh         OLF Trpc4     GSE131836 439 (raw downloaded) (CpG sites subsetted)
# 3: 5E_180925       Inh       MSN-D1 Hrh1     GSE134653 369 (raw downloaded) (CpG sites subsetted)
# 4: 4D_171219       Inh    MSN-D2 Slc24a2     GSE131406 353 (raw downloaded) (CpG sites subsetted)
# 5: 5J_190207       Inh    LSX-Inh Dock10     GSE135307 241 (raw downloaded) (CpG sites subsetted)
# 6: 4H_180911       Inh     PAL-Inh Ptprd     GSE134617 164 (raw downloaded) (CpG sites subsetted)
# 7: 5F_181220       Inh     MSN-D1 Plxnc1     GSE135169 158 (raw downloaded) (CpG sites subsetted)
# 8: 4F_180329       Inh D1L-Fstl4 Sipa1l2     GSE132445 117 (raw downloaded) (CpG sites subsetted)

as.list(cell_counts_inh_selected$GEO_accession[5]) %>% map(wrapper, mc.cores = 16) 
# Note: might cause memory exceeding with too many accessions processed simultaneously


##### particularly downloaded because of low coverage #####
# 4A_180205 GSE131766
# wrapper(geo_accession = "GSE131766", mc.cores = 16)
# wrapper(geo_accession = "GSE131354", mc.cores = 16)
# wrapper(geo_accession = "GSE131360", mc.cores = 16)
# wrapper(geo_accession = "GSE131427", mc.cores = 16)
# wrapper(geo_accession = "GSE131867", mc.cores = 16) 


##### other #####
# wrapper(geo_accession = "GSE135484", mc.cores = 16) 
# wrapper(geo_accession = "GSE158092", mc.cores = 16) 
# wrapper(geo_accession = "GSE135606", mc.cores = 16) 
# wrapper(geo_accession = "GSE135585", mc.cores = 16) 