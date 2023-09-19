source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")

# Individual cell files
outer_path <- "data/raw_counts/raw_counts_Luo2017_mice/"
cell_file_dirs <- list.files(outer_path)

# Take a subset of Luo 2017 for pilot study
metadata <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")
cell_counts <- metadata[, .(.N), by = .(Neuron_type1, Neuron_type3)] %>% arrange(desc(N)) 
#  3:   Excitatory          mL4 370

subtypes <- c("mL4")
sub_metadata <- metadata[Neuron_type3 %in% subtypes, .(sample, specie, Neuron_type1, Neuron_type2, Neuron_type3)]
nrow(sub_metadata) # = 370 cells in total

# cell file directoru
file_names <- paste0(sub_metadata$sample, "_indexed.tsv.gz")
stopifnot(all(file_names %in% cell_file_dirs))
sub_metadata$dir <- paste0(outer_path, file_names)

fwrite(sub_metadata, "data/interim/case_studies/luo2017mice_subset_hom/metadata_luo2017mice_subset_hom.csv")


