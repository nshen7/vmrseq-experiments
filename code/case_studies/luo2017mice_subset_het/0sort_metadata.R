.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
devtools::load_all("../vmrseq-package/vmrseq/")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

# Individual cell files
outer_path <- "data/raw_counts/raw_counts_Luo2021_mice/"
cell_file_dirs <- list.files(outer_path)

# Take a subset of Luo 2017 for pilot study
metadata <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")
cell_counts <- metadata[, .(.N), by = .(Neuron_type1, Neuron.type)] %>% arrange(desc(N)) 
# 5:   Excitatory       mDL-2 272
# 6:   Inhibitory         mPv 136
# 9:   Inhibitory      mSst-1 101
# 15:   Inhibitory      mSst-2  22
subtypes <- c("mDL-2", "mPv", "mSst-1", "mSst-2")
sub_metadata <- metadata[Neuron.type %in% subtypes, .(sample, specie, Neuron_type1, Neuron_type2, Neuron_type3, Neuron.type)]
nrow(sub_metadata) # = 531 cells in total

# cell file directoru
file_names <- paste0(sub_metadata$sample, "_indexed.tsv.gz")
stopifnot(all(file_names %in% cell_file_dirs))
sub_metadata$dir <- paste0(outer_path, file_names)

fwrite(sub_metadata, "data/interim/case_studies/luo2017mice_subset_het/metadata_luo2017mice_subset_het.csv")


