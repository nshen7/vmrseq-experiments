source("code/SETPATHS.R")
library(BiocParallel)
register(MulticoreParam(workers = 14))

source("code/case_studies/util_functions/4summarizeOutput.R")

read_dir <- "data/interim/case_studies/luo2017mice_full/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/"
if (!file.exists(write_dir)) dir.create(plot_dir)

# summarizeAllSites(read_dir, write_dir)
# summarizeOutputRegion(read_dir, write_dir, methods = c("vmrseq", "scbs", "smallwood", "scmet")) 
# summarizeOutputSite(read_dir, write_dir, methods = c("vmrseq", "scbs", "smallwood", "scmet")) 

# summarizeOutputRegion(read_dir, write_dir, methods = c("vmrseq"))

# summarizeOutputRegion(read_dir, write_dir, methods = c("scbs"))

# summarizeOutputRegion(read_dir, write_dir, methods = c("smallwood"))
summarizeOutputRegion(read_dir, write_dir, methods = c("smallwood_2kb"))

# summarizeOutputRegion(read_dir, write_dir, methods = c("scmet"))

