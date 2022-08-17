library(data.table)
setwd(here::here())

# ==== utils ====
file.size.head <- function(geo_accession, n){
  folder_dir <- paste0("data/raw_counts/raw_counts_Liu2021/",geo_accession,"_RAW/")
  # file_names0 <- list.files(folder_dir)
  file_names <- fread(paste0(folder_dir, "filelist.txt"))[-1,Name]
  file_dirs <- paste0(folder_dir, file_names)
  file_size <- file.size(file_dirs)
  file_size_head <- round(sort(file_size, decreasing = T)[1:n] / 2^20)
  return(file_size_head)
}

# file.mtime.head <- function(geo_accession, n){
#   folder_dir <- paste0("data/raw_counts/raw_counts_Liu2021/",geo_accession,"_RAW/")
#   # file_names0 <- list.files(folder_dir)
#   file_names <- fread(paste0(folder_dir, "filelist.txt"))[-1,Name]
#   file_dirs <- paste0(folder_dir, file_names)
#   file_mtime <- file.mtime(file_dirs)
#   file_mtime_head <- sort(file_mtime, decreasing = T)[1:n]
#   return(file_mtime_head)
# }


# ==== inspection ====

metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv")
GEO_all <- unique(metadata$GEO_accession); GEO_all <- GEO_all[!is.na(GEO_all)]

GEO_exist <- gsub("(.*)_RAW", "\\1", list.files("data/raw_counts/raw_counts_Liu2021/"))
GEO_exist <- GEO_exist[GEO_exist!="README.txt"]

## check if all GEO accessions have been downloaded
all(GEO_all %in% GEO_exist)
all(GEO_exist %in% GEO_all)

## check the maximum individual cell bedfile size is small (<= 60 MB),
## i.e., check if all GEO accessions have been fully subseted to only CpG sites
for (geo in GEO_exist) {
  cat(geo, file.size.head(geo_accession = geo, n = 1), "\n")
}

## check if bedfiles of all cells in metadata can be located
metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv") %>%
  filter(!is.na(SubType))
file_paths <- paste0("data/raw_counts/raw_counts_Liu2021/", 
                     metadata$GEO_accession, "_RAW/",
                     metadata$fullCEMBA)
file.exists(file_paths[1])
