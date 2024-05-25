source("code/SETPATHS.R")
source("code/case_studies/util_functions/2-4formatDataPca.R")

read_dir <- "data/interim/case_studies/luo2017mice_full/vmrseq/input"
write_dir <- "data/interim/case_studies/luo2017mice_full/pca/input"
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

formatData100kbins(read_dir, write_dir, chrs = paste0("chr", 1:19), bp_size = 100000, n_cores = 14)