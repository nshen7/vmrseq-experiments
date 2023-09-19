# ref: Smallwood, S., Lee, H., Angermueller, C. et al. https://doi.org/10.1038/nmeth.3035
source('code/SETPATHS.R')
library(parallel)
source("code/case_studies/util_functions/3-3runSmallwood.R")

read_dir <- "data/interim/case_studies/luo2017mice_subset_hom/smallwood/input/"
write_dir <- "data/interim/case_studies/luo2017mice_subset_hom/smallwood/output/"

runSmallwood(read_dir, write_dir, chrs = paste0("chr", 1:19), alpha = 0.05, n_cores = 14)