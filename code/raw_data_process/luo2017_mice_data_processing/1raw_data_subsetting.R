source("code/SETPATHS.R")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))

read_from <- "data/raw_counts/luo2017_from_geo/GSE97179"
write_to <- "data/raw_counts/raw_counts_Luo2017_mice/"

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

subsetAndStrandFlip_1cell <- function(write_to, cell_folder_name, chr_num, mc.cores) { 
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

# folder_list <- list.files(read_from)
# 
# ### folder_list_modified are samples names but is kept in the same order as folder_list
# folder_list_modified <- gsub("-", "_", folder_list)
# folder_list_modified <- gsub("(nuclei_.*)_.*_.*$", "\\1", folder_list_modified)
# folder_list_modified <- gsub("(.*_indexed).*", "\\1", folder_list_modified) 
# 
# metadata <- fread("data/metadata/metadata_luo2017/luo_supp_tables/NIHMS893063-supplement-Table_S1_csv.csv", skip = 1, header = T)
# all(metadata$Sample %in% folder_list_modified) # = TRUE
# 
# folder_list_mice <- folder_list[which(folder_list_modified %in% metadata$Sample)]
# 
# for (cell_folder_name in folder_list_mice[1348:length(folder_list_mice)]) {
#   subsetAndStrandFlip_1cell(cell_folder_name, chr_num = 1:19, mc.cores = 16)
# }

# ==== collapse counts of duplicated CpGs due to strand flipping ====

file_list <- list.files(write_to)
for (file in file_list) {
  # cell.df <- fread(here(write_to, file)) 
  # cell_collapsed.df <- cell.df %>%
  #   group_by(chr, strand, pos) %>%
  #   summarise(mc_count = sum(mc_count), 
  #             total = sum(total))
  # fwrite(cell_collapsed.df, here(write_to, file), sep = '\t')
  
  cell.df <- fread(here(write_to, file)) %>%
    mutate(strand = '*') %>%
    select(chr, pos, strand, mc_count, total)
  fwrite(cell.df, here(write_to, file), sep = '\t')
}


# ----- misc -----
# Distribution of proportion of sites with intermediate methylation level
countItmdProportion <- function(file) {
  cell.df <- fread(here(write_to, file))
  mf <- cell.df$mc_count/cell.df$total 
  return(sum(mf > 0 & mf < 1)/length(mf))
}
itmd_prop <- unlist(parallel::mclapply(list.files(write_to), countItmdProportion, mc.cores = 8))
quantile(itmd_prop)
#          0%         25%         50%         75%        100% 
# 0.001550225 0.003440104 0.004358966 0.005245920 0.021427392 
itmd_prop_all <- c(itmd_prop, itmd_prop_argelaguet) # concatenate with argelaguet study
quantile(itmd_prop_all)
#          0%         25%         50%         75%        100% 
# 0.000000000 0.003325643 0.004421386 0.005558797 0.322671732 


# Distribution of proportion of sites with <=2 reads
countUnder2readProportion <- function(file) {
  cell.df <- fread(here(write_to, file))
  return(sum(cell.df$total <= 2)/nrow(cell.df))
}
under2read_prop <- unlist(parallel::mclapply(list.files(write_to), countUnder2readProportion, mc.cores = 8))
quantile(under2read_prop)
#        0%       25%       50%       75%      100% 
# 0.7006976 0.9546422 0.9668191 0.9772215 0.9943490
under2read_prop_all <- c(under2read_prop, under2read_prop_argelaguet) # concatenate with argelaguet study
quantile(under2read_prop_all)
#        0%       25%       50%       75%      100% 
# 0.3783257 0.9413365 0.9620004 0.9754822 1.0000000 