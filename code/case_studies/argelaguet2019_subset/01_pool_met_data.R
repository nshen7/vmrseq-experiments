source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")

read_dir_met <- here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation_formatted', 'met_formatted')
write_dir <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', 'met', 'vmrseq', 'input')
if (!file.exists(write_dir)) dir.create(write_dir)

# Individual cell files
cell_file_names <- list.files(read_dir_met)

# Read in metadata
md <- fread(here('data', 'metadata', 'argelaguet2019', 'argelaguet2019_subset_met&rna_sample_metadata_processed.csv'))

# Sanity check whether all cells in metadata exist
stopifnot(all(na.omit(md$file_met) %in% cell_file_names))

## ---- Process methylation data ----
# Pool data
cell_file_dirs <- here(read_dir_met, md$file_met)
# data.pool(cellFiles = cell_file_dirs, writeDir = write_dir, chrNames = 1:2, colData = md)
# data.pool(cellFiles = cell_file_dirs, writeDir = write_dir, chrNames = 3:5, colData = md)
# data.pool(cellFiles = cell_file_dirs, writeDir = write_dir, chrNames = 6:9, colData = md)
# data.pool(cellFiles = cell_file_dirs, writeDir = write_dir, chrNames = 10:14, colData = md)
data.pool(cellFiles = cell_file_dirs, writeDir = write_dir, chrNames = 15:19, colData = md)
# 
