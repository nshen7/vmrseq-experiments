.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") 
source("code/case_studies/util_functions/2-2formatDataScmet.R")

read_dir <- "data/interim/case_studies/luo2017mice_subset_hom/vmrseq/input/"
write_dir <- "data/interim/case_studies/luo2017mice_subset_hom/scmet/input/"

formatDataScmet(read_dir, write_dir, chrs = paste0("chr", 1:19), bp_size = 20000, n_cores = 14)
