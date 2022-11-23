.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
devtools::load_all("../vmrseq-package/vmrseq/")

md <- fread("data/interim/case_studies/luo2017mice_subset_hom/metadata_luo2017mice_subset_hom.csv")
data.pool(cellFiles = md$dir,
          writeDir = "data/interim/case_studies/luo2017mice_subset_hom/vmrseq/input/",
          chrNames = paste0("chr", 1:19))
