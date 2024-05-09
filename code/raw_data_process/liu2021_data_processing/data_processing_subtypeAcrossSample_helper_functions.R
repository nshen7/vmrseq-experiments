suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(parallel))

source("code/SETPATHS.R")

metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv") %>%
  filter(!is.na(GEO_accession) & !is.na(SubType) & !is.na(FilePath))

# Estimate a maximum coordinate for each chr
sc.exp <- fread("data/raw_counts/raw_counts_Liu2021/GSE135269_RAW/GSM4001537_allc_CEMBA190205-6D-1-CEMBA190205-6D-2-E3_ad007.tsv.gz")
chr_pos <- sc.exp %>% 
  group_by(chr) %>% 
  summarise(max_pos = max(pos)) %>% 
  mutate(max_pos = max_pos + 1e6) %>%
  as.data.table()

summarize_1chr <- function(chr_name, file_dirs) {
  
  chr_len <- chr_pos[chr == chr_name, 'max_pos'] %>% unlist %>% unname
  smr.dt <- data.table(cell_meth = rep(0, chr_len),
                       cell_cov = rep(0, chr_len))
  
  ## Sequentially process each cell
  for (file_dir in file_dirs) {
    chr.dt <- fread(file_dir) %>% 
      filter(chr == chr_name) %>% 
      mutate(bool_meth = round(meth / total))
    
    stopifnot("Non-0/1 methylation value." = all(chr.dt$bool_meth %in% 0:1))
    stopifnot("Maximum genomic position exceeded." = max(chr.dt$pos) < chr_len)
    
    smr.dt[chr.dt$pos, 'cell_cov'] <- smr.dt[chr.dt$pos, 'cell_cov'] + 1
    smr.dt[chr.dt$pos, 'cell_meth'] <- smr.dt[chr.dt$pos, 'cell_meth'] + chr.dt$bool_meth
  }
  
  index <- which(smr.dt$cell_cov > 0)
  smr.dt <- smr.dt %>% 
    dplyr::slice(index) %>% 
    dplyr::mutate(chr = chr_name, pos = index) %>% 
    dplyr::select(chr, pos, cell_meth, cell_cov)
  
  return(smr.dt)
}

wrapper <- function(subtype, chr_names = paste0("chr", 1:19), mc.cores) {
  print(subtype)
  file_dirs <- with(metadata,  FilePath[which(SubType == subtype)])
  file_dirs <- file_dirs[!is.na(file_dirs)]
  
  dt <- do.call(rbind, mclapply(chr_names, 
                                function(.x) summarize_1chr(.x, file_dirs),
                                mc.cores = mc.cores))
  
  if (all(chr_names == paste0("chr", 1:19))) {
    write_dir <- paste0("data/processed/summarized_liu2021/summarized_subtype_", sub(" ", "_", subtype), "_", length(file_dirs), "cells.csv")
  } else {
    write_dir <- paste0("data/processed/summarized_liu2021/summarized_subtype_", sub(" ", "_", subtype), "_", length(file_dirs), "cells_", chr_names[1],"_to_",chr_names[length(chr_names)],".csv")
  }
  
  fwrite(dt, write_dir, quote = F)
  R.utils::gzip(write_dir, remove = T, overwrite = T)
  
}

