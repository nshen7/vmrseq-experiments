library(data.table)
library(stringr)

metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv")
GEO_exist <- gsub("(.*)_RAW", "\\1", list.files("data/raw_counts/raw_counts_Liu2021/")) 
GEO_all <- unique(metadata$GEO_accessio); GEO_all <- GEO_all[!is.na(GEO_all)]

GEO_tbd <- GEO_all[!GEO_all %in% GEO_exist]

script_dir <- "code/liu2021_data_processing/misc/download_GEO_bash_script.txt"
if(file.exists(script_dir)) unlink(script_dir)
for (GEO in GEO_tbd) {
  folder <- paste0(substr(GEO, 1, 6), "nnn")
  print(folder)
  cmd <- paste0("### Downloaded? 
cd /scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/data/raw_counts/raw_counts_Liu2021
n='",GEO,"'
wget \"https://ftp.ncbi.nih.gov/geo/series/",folder,"/$n/suppl/${n}_RAW.tar\"
mkdir \"${n}_RAW\"
tar -xvf \"${n}_RAW.tar\" -C \"${n}_RAW\"
wget \"https://ftp.ncbi.nih.gov/geo/series/",folder,"/$n/suppl/filelist.txt\" -P \"${n}_RAW\"
rm \"${n}_RAW.tar\"

")
  write(cmd, script_dir, append = T)
}
