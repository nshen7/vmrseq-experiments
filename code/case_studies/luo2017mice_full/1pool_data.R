.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
devtools::load_all("../vmrseq-package/vmrseq/")

# Individual cell files
outer_path <- "data/raw_counts/raw_counts_Luo2021_mice/"
cell_file_dirs <- list.files(outer_path)

# Sanity check whether all cells in metadata exist
md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")
file_names <- paste0(md$sample, "_indexed.tsv.gz")
stopifnot(all(file_names %in% cell_file_dirs))
md$dir <- paste0(outer_path, file_names)

# Pool data
write_dir <- "data/interim/case_studies/luo2017mice_full/vmrseq/input/"
if (!file.exists(write_dir)) dir.create(write_dir)

# data.pool(cellFiles = md$dir, writeDir = write_dir, chrNames = paste0("chr", 1:19))
# data.pool(cellFiles = md$dir, writeDir = write_dir, chrNames = paste0("chr", 10:12))
# data.pool(cellFiles = md$dir, writeDir = write_dir, chrNames = paste0("chr", 13:15))
data.pool(cellFiles = md$dir, writeDir = write_dir, chrNames = paste0("chr", 16:19))
