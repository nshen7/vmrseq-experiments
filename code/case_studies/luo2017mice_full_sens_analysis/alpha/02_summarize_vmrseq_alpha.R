source("code/SETPATHS.R")
library(BiocParallel)
register(MulticoreParam(workers = 14))

source("code/case_studies/util_functions/4summarizeOutput.R")

read_dir_1 <- "data/interim/case_studies/luo2017mice_full/"
read_dir_2 <- "data/interim/case_studies/luo2017mice_full_sens_analysis/alpha"
write_dir <- "data/interim/case_studies/luo2017mice_full_sens_analysis/result_summary/alpha"
if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- "plots/case_studies/luo2017mice_full_sens_analysis"
if (!file.exists(write_dir)) dir.create(plot_dir)


# ---- utils ----
### Summarize region methylation
summarizeOutputRegion <- function(read_dir_1, read_dir_2, write_dir, methods, chr_names = paste0('chr', 1:19)) {
  
  if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)
  
  se_dirs <- here(read_dir_1, "vmrseq", "input", chr_names)
  
  ### vmrseq
  if ("vmrseq" %in% methods) {
    dirs_vseq <- here(read_dir_2, "vmrseq", "output", paste0(chr_names, ".rds"))
    
    # Concatenate results from 19 chromosomes
    res_vseq <- list(
      vmr.ranges = do.call(c, lapply(dirs_vseq, function(dir) if (file.exists(dir)) readRDS(dir)$vmr.ranges else NULL)),
      cr.ranges = do.call(c, lapply(dirs_vseq, function(dir) if (file.exists(dir)) readRDS(dir)$cr.ranges else NULL))
    )
    rs_vseq1.se <- getRegionSummary(
      se_dirs = se_dirs,
      regions.gr = res_vseq$vmr.ranges
    )
    saveHDF5SummarizedExperiment(rs_vseq1.se, here(write_dir, "vmrseq_regionSummary_vmrs"), replace=TRUE)
    rs_vseq2.se <- getRegionSummary(
      se_dirs = se_dirs, 
      regions.gr = res_vseq$cr.ranges
    )
    saveHDF5SummarizedExperiment(rs_vseq2.se, here(write_dir, "vmrseq_regionSummary_crs"), replace=TRUE)
  }
}

# ---- run ----
alpha <- 0.01
# alpha <- 0.025
# alpha <- 0.05
# alpha <- 0.10
summarizeOutputRegion(read_dir_1 = read_dir_1, 
                      read_dir_2 = here(read_dir_2, paste0('alpha', alpha)), 
                      write_dir = here(write_dir, paste0('alpha', alpha)), 
                      methods = c("vmrseq"))
