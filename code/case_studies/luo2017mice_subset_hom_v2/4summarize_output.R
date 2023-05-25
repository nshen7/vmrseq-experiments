source("code/SETPATHS.R")
library(BiocParallel)
register(MulticoreParam(workers = 8))

source("code/case_studies/util_functions/4summarizeOutput.R")

read_dir <- "data/interim/case_studies/luo2017mice_subset_hom_v2/"
write_dir <- "data/interim/case_studies/luo2017mice_subset_hom_v2/result_summary/"
if (!file.exists(write_dir)) dir.create(write_dir)

summarizeOutputRegion(read_dir, write_dir, methods = c("vmrseq"))
summarizeOutputSite(read_dir, write_dir, methods = c("vmrseq"))

