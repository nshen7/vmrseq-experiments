source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")

read_dir <- here('data', 'interim', 'case_studies', 'luo2017mice_subset_hom_v2')
write_dir <- here('data', 'interim', 'case_studies', 'luo2017mice_subset_hom_v2', 'vmrseq')
if (!file.exists(write_dir)) dir.create(write_dir)

md <- fread(here(read_dir, "metadata_luo2017mice_subset_hom.csv"))
data.pool(cellFiles = md$dir,
          writeDir = here(write_dir, "input"),
          chrNames = paste0("chr", 1:19))
