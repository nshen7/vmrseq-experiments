.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(DelayedArray))
suppressPackageStartupMessages(library(DelayedMatrixStats))
suppressPackageStartupMessages(library(SummarizedExperiment))


# Use real data (subtype 'IT-L23 Cux1' in Liu2021 data) as null
data_source <- "liu2021"; subtype <- "IT-L23 Cux1"
chromosome <- "chr1"
folder <- paste0("data/interim/sim_studies/real/", 
                 data_source, "_", sub(" ", "_", subtype), "_", chromosome)
if (!dir.exists(folder)) dir.create(folder)

# Get file directories for inidividual cells
metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv") %>%
  filter(!is.na(GEO_accession) & !is.na(SubType) & !is.na(FilePath))
file_dirs <- with(metadata,  FilePath[which(SubType == subtype)])
rm(metadata)

# Function for process 1 cell
processCell <- function(i, select_chrs) {
  file_dir <- file_dirs[i]
  cell_df <- fread(file_dir) %>%
    filter(chr %in% select_chrs) %>%
    filter(meth/total == 1 | meth/total == 0) %>%
    filter(!duplicated(pos)) %>%
    mutate(bool = round(meth/total)) %>%
    dplyr::select(chr, pos, bool)
  # colnames(cell_df)[3] <- paste0("cell_",i)
  return(cell_df)
}


# Get all possible genomic coordinates in chromosome
pos0 <- processCell(1, chromosome)$pos
for (i in 2:length(file_dirs)) {
  pos_temp <- processCell(i, chromosome)$pos
  pos_new <- pos_temp[which(!pos_temp %in% pos0)]
  pos0 <- c(pos0, pos_new)
  cat(i, " ")
}
pos0 <- sort(pos0)

pos_dir <- paste0(folder, "/pos0.txt")
fwrite(data.frame(pos0), pos_dir, col.names = F)
R.utils::gzip(pos_dir, remove = T, overwrite = T)

# pos0 <- fread(paste0(folder, "/pos0.txt.gz")) %>% unlist %>% unname


# Store individual cell files with NAs added in to uncovered sites
fillCell <- function(i) {
  # Read in cell info
  cell_df <- processCell(i, chromosome)
  
  # Fill in NA for uncovered position
  filled <- rep(NA, length(pos0))
  ind <- which(pos0 %in% cell_df$pos)
  if(length(ind) != nrow(cell_df)) stop("'pos0' does not contain all possible positions")
  filled[ind] <- cell_df$bool
  
  # Write out filled info
  write_dir <- paste0(folder, "/cell_", i, ".txt")
  fwrite(as.data.frame(filled), write_dir, col.names = F, na = "NA", quote = F)
  R.utils::gzip(write_dir, remove = T, overwrite = T)
}

bplapply(1:length(file_dirs), fillCell)
