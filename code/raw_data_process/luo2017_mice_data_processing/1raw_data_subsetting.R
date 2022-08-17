suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")

read_from <- "../../DXM_extend_chr1/data/raw_counts/"
write_to <- "data/raw_counts/raw_counts_Luo2021_mice/"

# ==== helper functions ====
subsetAndStrandFlip_1chr <- function(chr_file) {
  chr <- fread(chr_file, header = T,
               colClasses = c("character", "integer", "character", "character","integer", "integer", "integer"))
  # colnames inherited from data: chr, pos, strand, mc_class, mc_count, total, methylated
  chr$methylated <- NULL
  
  ## remove all non-CpG calls.
  chr <- chr[grep(pattern="CG.", x=mc_class), ]
  chr$mc_class <- NULL
  
  chr$chr <- paste("chr", chr$chr, sep = "")
  return(chr)
}

subsetAndStrandFlip_1cell <- function(cell_folder_name, chr_num, mc.cores) { 
  # `cell_folder_name` -- folder for chromosome files for 1 cell
  # `chr_num` -- chromosome numbers, vector of integer number and/or single characters (e.g., c(1,2,3,X,Y))
  
  file_names <- list.files(paste0(read_from, cell_folder_name))
  selected_files <- file_names[gsub(".*_(.*).tsv.gz$", "\\1",file_names) %in% as.character(chr_num)]
  chr_files <- paste0(read_from, cell_folder_name, "/", selected_files) # chromosome file path
  
  cell <- do.call(rbind, lapply(chr_files, subsetAndStrandFlip_1chr))
  cell$pos[cell$strand == "-"] <- cell$pos[cell$strand == "-"] - 1
  cell$strand <- "*"
  cat("End processing cell", cell_folder_name, "\n")
  
  write_dir <- paste0(write_to, cell_folder_name, ".tsv")
  write_tsv(cell, write_dir)
  R.utils::gzip(write_dir, remove = T, overwrite = T)
  cat("End writing file for cell", cell_folder_name, "\n")
}


# ==== raw data subsetting ====

folder_list <- list.files(read_from)

### folder_list_modified are samples names but is kept in the same order as folder_list
folder_list_modified <- gsub("-", "_", folder_list)
folder_list_modified <- gsub("(nuclei_.*)_.*_.*$", "\\1", folder_list_modified)
folder_list_modified <- gsub("(.*_indexed).*", "\\1", folder_list_modified) 

metadata <- fread("data/metadata/metadata_luo2017/NIHMS893063-supplement-Table_S1_csv.csv", skip = 1, header = T)
all(metadata$Sample %in% folder_list_modified) # = TRUE

folder_list_mice <- folder_list[which(folder_list_modified %in% metadata$Sample)]

for (cell_folder_name in folder_list_mice[1348:length(folder_list_mice)]) {
  subsetAndStrandFlip_1cell(cell_folder_name, chr_num = 1:19, mc.cores = 16)
}
