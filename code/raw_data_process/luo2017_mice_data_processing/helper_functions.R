suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(bsseq))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(DelayedArray))

# ==== main function ====

### convert raw count data to SummerizedExperiment object
rawCountToSE <- function(cell_file_dirs, metadata, cell_ids, write_dir, mc.cores){
  
  if(length(cell_file_dirs) != length(cell_ids)) stop("Number of file names not equal to number of sample names.")
  if(length(cell_file_dirs) != nrow(metadata)) stop("Number of file names not equal to number of rows in metadata.")
  
  dat.bs <- bsseq::combineList(mclapply(1:length(cell_ids),
                                        function(i) with(fread(cell_file_dirs[i]),
                                                         BSseq(gr = GRanges(seqnames = Rle(chr),
                                                                            ranges = IRanges(start = pos, end = pos),
                                                                            strand = strand),
                                                                            M = matrix(mc_count),
                                                                            Cov = matrix(total),
                                                                            sampleNames = cell_ids[i]))
                                        ,
                                        mc.cores = mc.cores
                                        )
                               )
  print("Finished combining cells.")
  
  dat.se <- SummarizedExperiment(assays = list("Cov" = getCoverage(dat.bs), "M" = getCoverage(dat.bs, type = "M")),
                                 rowRanges = granges(dat.bs),
                                 colData = metadata)
  print("Finished creating SummarizedExperiment object.")
  saveHDF5SummarizedExperiment(dat.se, dir = write_dir, replace = T)
  print("Finished writing HDF5 file.")
}


# ==== sort out sample types from supp table S2 of scMET paper ====

wrapper <- function(subtype, mc.cores){

  metadata <- fread("../../DXM_extend_chr1/data/metadata/sample_info_processed.csv")
  
  ### input arguments
  outer_path <- "data/raw_counts/raw_counts_Luo2021_mice/"
  cell_file_dirs <- list.files(outer_path)
  sub_metadata <- metadata[Neuron_type3 == subtype, .(sample, specie, Neuron_type1, Neuron_type2, Neuron_type3)]
  
  stopifnot(all(paste0(sub_metadata$sample, "_indexed.tsv.gz") %in% cell_file_dirs))
  sub_metadata$folder_dir <- paste0(outer_path, sub_metadata$sample, "_indexed.tsv.gz")
  
  n.cells <- nrow(sub_metadata)
  cat("Total", n.cells, "cells are in subtype", subtype, "\n")

  ### write out bsseq objects 
  rawCountToSE(cell_file_dirs = sub_metadata$folder_dir, 
               metadata = sub_metadata,  
               cell_ids = sub_metadata$sample,
               write_dir = paste0("data/processed/processed_luo2017_mice/subtype_", sub("/", "", subtype), "_", n.cells, "cells"),
               mc.cores = mc.cores)
}
