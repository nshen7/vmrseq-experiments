source("code/SETPATHS.R")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(SummarizedExperiment))

write_dir <- "../vmrseq-package/vmrseq-workflow-vignette/data/example/"
md <- fread("data/interim/case_studies/luo2017mice_subset_het/metadata_luo2017mice_subset_het.csv") %>%
  dplyr::select(sample, specie, Neuron_type1, Neuron.type) %>%
  dplyr::rename(broad_class = Neuron_type1, subtype = Neuron.type)

chr <- "chr1"
SE <- loadHDF5SummarizedExperiment(paste0("data/interim/case_studies/luo2017mice_subset_het/vmrseq/input/", chr))
colData(SE) <- cbind(colData(SE), md)
saveHDF5SummarizedExperiment(SE, here(write_dir, chr), replace=TRUE)


