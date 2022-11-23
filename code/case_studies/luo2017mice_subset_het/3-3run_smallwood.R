# ref: Smallwood, S., Lee, H., Angermueller, C. et al. https://doi.org/10.1038/nmeth.3035
.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
source("code/case_studies/util_functions/3-3runSmallwood.R")

read_dir <- "data/interim/case_studies/luo2017mice_subset_het/smallwood/input/"
write_dir <- "data/interim/case_studies/luo2017mice_subset_het/smallwood/output/"

runSmallwood(read_dir, write_dir, chrs = paste0("chr", 1:19), alpha = 0.05, n_cores = 14)