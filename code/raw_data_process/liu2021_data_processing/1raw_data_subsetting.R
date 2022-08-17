## !!! This script processes raw bed files downloaded from GEO and extract out the CpGs and remove all non-CpGs.

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")

# ==== helper functions ====
subsetAndStrandFlip <- function(file_dir) {
  if (file.size(file_dir) < 30*2^20) return()
  # if (ncol(fread(file_dir)) == 6) return()
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
  folder_dir <- paste0("data/raw_counts/raw_counts_Liu2021/",geo_accession,"_RAW/")
  file_names <- fread(paste0(folder_dir, "filelist.txt"))[-1,Name]
  file_dirs <- paste0(folder_dir, file_names)
  mclapply(file_dirs, subsetAndStrandFlip, mc.cores = mc.cores)
  # mclapply(file_dirs[1], subsetAndStrandFlip, mc.cores = mc.cores) ### for testing purpose
  cat("End processing", geo_accession, "\n")
}



# ==== raw data subsetting ====
# wrapper(geo_accession = "GSE131519", mc.cores = 16)
# wrapper(geo_accession = "GSE131557", mc.cores = 16)
# wrapper(geo_accession = "GSE131951", mc.cores = 16)
# wrapper(geo_accession = "GSE131554", mc.cores = 16)
# wrapper(geo_accession = "GSE132060", mc.cores = 16)

# wrapper(geo_accession = "GSE132114", mc.cores = 16)
# wrapper(geo_accession = "GSE132126", mc.cores = 16)
# wrapper(geo_accession = "GSE132498", mc.cores = 16)
# wrapper(geo_accession = "GSE132503", mc.cores = 16)
# wrapper(geo_accession = "GSE132463", mc.cores = 16)

# wrapper(geo_accession = "GSE131913", mc.cores = 16)
# wrapper(geo_accession = "GSE132101", mc.cores = 16)
# wrapper(geo_accession = "GSE132639", mc.cores = 16)

# wrapper(geo_accession = "GSE131757", mc.cores = 16)
# wrapper(geo_accession = "GSE131570", mc.cores = 16)
# wrapper(geo_accession = "GSE132541", mc.cores = 16)

# wrapper(geo_accession = "GSE131192", mc.cores = 16) 
# wrapper(geo_accession = "GSE131326", mc.cores = 16) 
# wrapper(geo_accession = "GSE132571", mc.cores = 16)
# wrapper(geo_accession = "GSE132135", mc.cores = 16) 

# wrapper(geo_accession = "GSE131771", mc.cores = 16)
# wrapper(geo_accession = "GSE132484", mc.cores = 16) 
# wrapper(geo_accession = "GSE132405", mc.cores = 16)

# wrapper(geo_accession = "GSE131430", mc.cores = 16)
# wrapper(geo_accession = "GSE132601", mc.cores = 16)
# wrapper(geo_accession = "GSE132609", mc.cores = 16)

# wrapper(geo_accession = "GSE131509", mc.cores = 16)
# wrapper(geo_accession = "GSE131393", mc.cores = 16) 
# wrapper(geo_accession = "GSE131501", mc.cores = 16) 

# wrapper(geo_accession = "GSE132456", mc.cores = 16)
# wrapper(geo_accession = "GSE134635", mc.cores = 16)
# wrapper(geo_accession = "GSE134940", mc.cores = 16)
# wrapper(geo_accession = "GSE134490", mc.cores = 16)

# wrapper(geo_accession = "GSE135260", mc.cores = 16)
# wrapper(geo_accession = "GSE135641", mc.cores = 16)

# wrapper(geo_accession = "GSE132676", mc.cores = 16)
# wrapper(geo_accession = "GSE134750", mc.cores = 16)
# wrapper(geo_accession = "GSE134967", mc.cores = 16)

# wrapper(geo_accession = "GSE132700", mc.cores = 16)
# wrapper(geo_accession = "GSE132710", mc.cores = 16)
# wrapper(geo_accession = "GSE134958", mc.cores = 16)
# wrapper(geo_accession = "GSE135139", mc.cores = 16)

# wrapper(geo_accession = "GSE134793", mc.cores = 16)
# wrapper(geo_accession = "GSE134507", mc.cores = 16)
# wrapper(geo_accession = "GSE134593", mc.cores = 16) 
# wrapper(geo_accession = "GSE135203", mc.cores = 16)

# wrapper(geo_accession = "GSE134882", mc.cores = 16)

# wrapper(geo_accession = "GSE134899", mc.cores = 16)
# wrapper(geo_accession = "GSE135234", mc.cores = 16)
# wrapper(geo_accession = "GSE135248", mc.cores = 16)
# wrapper(geo_accession = "GSE135252", mc.cores = 16)

# wrapper(geo_accession = "GSE135269", mc.cores = 16)
# wrapper(geo_accession = "GSE132621", mc.cores = 16)
# wrapper(geo_accession = "GSE132629", mc.cores = 16)
# wrapper(geo_accession = "GSE134474", mc.cores = 16)

# wrapper(geo_accession = "GSE134484", mc.cores = 16)
wrapper(geo_accession = "GSE158043", mc.cores = 16)
wrapper(geo_accession = "GSE158154", mc.cores = 16)



