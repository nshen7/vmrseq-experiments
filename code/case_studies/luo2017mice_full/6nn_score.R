.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
source("code/case_studies/util_functions/6nnScore.R")

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")

## Regional mean methyl
nnScoreMethod <- function(method, top_n_regions = NULL, k = 100, theta = 0.9) {
  
  name_seg <- ifelse(is.null(top_n_regions), yes = "", no = paste0("_top", top_n_regions, "regions"))
  d_mat <- fread(
    paste0(write_dir, "dissimilarity_matrix_regional_methyl_", method, name_seg, "_seed2010.txt.gz"), 
    drop = 1
  )
  
  stopifnot(all(md$sample == names(d_mat)))
  true_clust <- md$Neuron.type
  
  return(nnScore(d_mat, true_clust, k, theta))
}

nnScoreMethod(method = "vseq") # = 0.5236233
nnScoreMethod(method = "vseq_cr") # = 0.4991854
nnScoreMethod(method = "scbs") # = 0.4095797
nnScoreMethod(method = "smwd") # =  0
# nnScoreMethod(method = "scmet") # = 

nnScoreMethod(method = "vseq", top_n_regions = 1000) # = 0.4301075
nnScoreMethod(method = "scbs", top_n_regions = 1000) # = 0.3515803
nnScoreMethod(method = "smwd", top_n_regions = 1000) # = 0
# nnScoreMethod(method = "scmet", top_n_regions = 1000) # = 

nnScoreMethod(method = "vseq", top_n_regions = 300) # = 0.3796025
nnScoreMethod(method = "scbs", top_n_regions = 300) # = 0.230694
nnScoreMethod(method = "smwd", top_n_regions = 300) # = 0
# nnScoreMethod(method = "scmet", top_n_regions = 300) # = 

