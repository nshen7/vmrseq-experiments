.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
source("code/case_studies/util_functions/6nnScore.R")

read_dir <- "data/interim/case_studies/luo2017mice_subset_het/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_subset_het/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_subset_het/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

## Regional mean methyl
nnScoreMethod <- function(method, top_n_regions = NULL, k = 20, theta = 0.9) {
  
  name_seg <- ifelse(is.null(top_n_regions), yes = "", no = paste0("_top", top_n_regions, "regions"))
  d_mat <- fread(
    paste0(write_dir, "dissimilarity_matrix_regional_methyl_", method, name_seg, "_seed2010.txt.gz"), 
    drop = 1
  )
  
  md <- fread("data/interim/case_studies/luo2017mice_subset_het/metadata_luo2017mice_subset_het.csv")
  stopifnot(all(md$sample == names(d_mat)))
  true_clust <- md$Neuron.type
  
  return(nnScore(d_mat, true_clust, k, theta))
}

nnScoreMethod(method = "vseq") # = 0.8436911
nnScoreMethod(method = "vseq_cr") # = 0.8248588
nnScoreMethod(method = "scbs") # = 0.7815443
nnScoreMethod(method = "smwd") # = 0.5913371
nnScoreMethod(method = "scmet") # = 0

nnScoreMethod(method = "vseq", top_n_regions = 1000) # = 0.7853107 
nnScoreMethod(method = "scbs", top_n_regions = 1000) # = 0.6290019
nnScoreMethod(method = "smwd", top_n_regions = 1000) # = 0.3352166
nnScoreMethod(method = "scmet", top_n_regions = 1000) # = 0.06967985

nnScoreMethod(method = "vseq", top_n_regions = 300) # = 0.6741996
nnScoreMethod(method = "scbs", top_n_regions = 300) # = 0.5235405
nnScoreMethod(method = "smwd", top_n_regions = 300) # = 0.1280603
nnScoreMethod(method = "scmet", top_n_regions = 300) # = 0.04143126

