.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
library(BiocParallel)
register(MulticoreParam(workers = 14))

source("code/case_studies/util_functions/4summarizeOutput.R")

read_dir <- "data/interim/case_studies/luo2017mice_subset_het/"
write_dir <- "data/interim/case_studies/luo2017mice_subset_het/result_summary/"
if (!file.exists(write_dir)) dir.create(write_dir)

# summarizeAllSites(read_dir, write_dir)
# summarizeOutputRegion(read_dir, write_dir, methods = c("vmrseq", "scbs", "smallwood")) 
# summarizeOutputSite(read_dir, write_dir, methods = c("vmrseq", "scbs", "smallwood")) 

summarizeOutputSite(read_dir, write_dir, methods = c("scmet"))
summarizeOutputRegion(read_dir, write_dir, methods = c("scmet"))


