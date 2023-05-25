source("code/SETPATHS.R")
source("code/case_studies/util_functions/2-1formatDataScbs.R")

write_dir <- "data/interim/case_studies/luo2017mice_full/scbs/input/"
if (!file.exists(write_dir)) dir.create(write_dir)

read_dir <- "data/interim/case_studies/luo2017mice_full/vmrseq/input/"
formatDataScbs(read_dir, write_dir, chrs = paste0("chr", 1:19), n_cores = 14)
