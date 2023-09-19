source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")

md <- fread("data/interim/case_studies/luo2017mice_subset_hom/metadata_luo2017mice_subset_hom.csv")
data.pool(cellFiles = md$dir,
          sep = '\t',
          writeDir = "data/interim/case_studies/luo2017mice_subset_hom/vmrseq/input/",
          chrNames = paste0("chr", 1:19))
