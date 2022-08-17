suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(bsseq))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(DelayedArray))

# ==== sub function ====
convertToBSseq <- function(file_dir, cell_id) {
  cell <- fread(
    file_dir,
    header = T,
    select = 1:6,
    colClasses = c("character", "integer", "character", "character",
                   "integer", "integer"))
  
  cell.bs <- with(cell, 
                  BSseq(gr = GRanges(seqnames = Rle(chr), 
                                     ranges = IRanges(start = pos, end = pos),
                                     strand = strand),
                        M = matrix(meth),
                        Cov = matrix(total),
                        sampleNames = cell_id))
  rm(cell)
  return(cell.bs)
}

# ==== main function ====
### convert raw count data to SummarizedExperiment object
rawCountToSE <- function(file_dirs, metadata, cell_ids, write_dir, mc.cores){
  
  if(length(file_dirs) != length(cell_ids)) stop("Number of file names not equal to number of sample names.")
  if(length(file_dirs) != nrow(metadata)) stop("Number of file names not equal to number of rows in metadata.")
  
  dat.bs = bsseq::combineList(mclapply(1:length(cell_ids),
                                       function(i) convertToBSseq(file_dirs[i], cell_ids[i]),
                                       mc.cores = mc.cores))
  # rm(dat_list.bs)
  print("Finished combining cells.")
  
  dat.se <- SummarizedExperiment(assays = list("Cov" = getCoverage(dat.bs), "M" = getCoverage(dat.bs, type = "M")),
                                 rowRanges = granges(dat.bs),
                                 colData = metadata)
  print("Finished creating SummarizedExperiment object.")
  saveHDF5SummarizedExperiment(dat.se, dir = write_dir, replace = T)
  print("Finished writing HDF5 file.")
}


# ==== wrapper function for processing and writing files ====
gsub2 <- function(x, ...) gsub(..., x = x)

wrapper <- function(sample_name, geo_accession, subtype, mc.cores){
  
  ### input arguments
  folder_dir <- paste0("data/raw_counts/raw_counts_Liu2021/", geo_accession, "_RAW/")
  file_names <- fread(paste0(folder_dir, "filelist.txt"))[-1,Name]
  
  #### modify file name format for matching metadata entries
  file_subCEMBA <- file_names %>% 
    gsub2(pattern = "ad", replacement = "AD") %>%
    gsub2(pattern = ".*(CEMBA[^_].*CEMBA[^_]..*AD\\d{3}).*", replacement = "\\1") %>%
    gsub2(pattern = "-", replacement = "_")
  
  ### check if matching succeeded
  metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv")
  print("Checking if file names are successfully matched to metadata entries.")
  n_matched <- sum(file_subCEMBA %in% metadata$subCEMBA)
  cat(n_matched,"cells out of", length(file_names), "from sample", sample_name,"are matched to metadata.\n")
  
  ### change to another matching rule if failed (to accommodate samples such as: sample_name = "4H_180911", geo_accession = "GSE134617")
  if(n_matched == 0){
    print("First attempt of matching failed. Trying alternative matching rule.")
    file_subCEMBA <- file_names %>% 
      gsub2(pattern = "ad", replacement = "AD") %>%
      gsub2(pattern = ".*(allc_.*AD\\d{3})", replacement = "\\1") %>%
      gsub2(pattern = "-", replacement = "_")
    print("Checking again if file names are successfully matched to metadata entries.")
    n_matched <- sum(file_subCEMBA %in% metadata$subCEMBA)
    cat(n_matched,"cells out of", length(file_names), "are matched to metadata.\n")
  }
  
  dt_files <- data.table(file_names, file_subCEMBA)
  if(n_matched < length(file_names)) dt_files <- dt_files[file_subCEMBA %in% metadata$subCEMBA]

  sub_metadata <- merge.data.table(dt_files, metadata,
                                   by.x = "file_subCEMBA", by.y = "subCEMBA",
                                   all.x = F, all.y = F)[SubType == subtype & Sample == sample_name & `Pass QC` == TRUE,
                                                         .(file_names, CellID, Sample, RegionName, MajorRegion, SubRegion, CellClass, MajorType, SubType)]
  
  n_cells <- nrow(sub_metadata)
  cat("Total", n_cells, "cells are in sample", sample_name, "subtype", subtype, "\n")
  
  ### write out bsseq objects
  file_dirs <- paste0(folder_dir, sub_metadata$file_names)
  rawCountToSE(file_dirs = file_dirs, metadata = sub_metadata,  cell_ids = sub_metadata$CellID,
               write_dir = paste0("data/processed/processed_liu2021/sample",sample_name,
                                  "_",geo_accession, "_subtype_", sub(" ", "_", subtype), "_", n_cells, "cells"),
               mc.cores = mc.cores)
  print("\n")
}
