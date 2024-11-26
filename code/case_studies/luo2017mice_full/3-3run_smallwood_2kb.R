# ref: Smallwood, S., Lee, H., Angermueller, C. et al. https://doi.org/10.1038/nmeth.3035
source("code/SETPATHS.R")
source("code/case_studies/util_functions/3-3runSmallwood.R")

read_dir <- "data/interim/case_studies/luo2017mice_full/smallwood_2kb/input/"
write_dir <- "data/interim/case_studies/luo2017mice_full/smallwood_2kb/output/"
if (!file.exists(write_dir)) dir.create(write_dir, recursive = TRUE)

runSmallwood(read_dir, write_dir, chrs = paste0("chr", 1:19), alpha = 0.05, n_cores = 22)