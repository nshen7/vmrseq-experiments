source("code/SETPATHS.R")
source("code/case_studies/util_functions/2-2formatDataScmet.R")

read_dir <- "data/interim/case_studies/luo2017mice_full/vmrseq/input/"
write_dir <- "data/interim/case_studies/luo2017mice_full/scmet/input/"
if (!file.exists(write_dir)) dir.create(write_dir)

formatDataScmet(read_dir, write_dir, chrs = paste0("chr", 1:19), bp_size = 20000, n_cores = 14)