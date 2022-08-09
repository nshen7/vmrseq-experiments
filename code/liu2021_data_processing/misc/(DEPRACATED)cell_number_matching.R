library(data.table)
library(parallel)

setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
file_list <- fread("data/raw_counts/raw_counts_Liu2021/GSE130553_RAW/filelist.txt")[-1,]

# ==== import metadata ====
md <- as.data.table(readxl::read_xlsx("data/metadata/metadata_liu2021/Liu2021_cell_metadata.xlsx", skip = 15))
colnames(md)[1] <- "CellName"

## PS: GSE130553 corresponds to sample CEMBA180226_1A
md_1A226 <- md[Sample=="1A_180226"]

md_1A226 <- data.table(cellNumber = as.integer(gsub(pattern = "1A_M_(.*)", replacement = "\\1", x = md_1A226$CellName)), md_1A226)

# ==== summary data of first 100 files in GSE130553 ====
calFracs <- function(file_dir){
  dt <- fread(file_dir, select = c(1,4,5,6))
  colnames(dt) <- c("chr","context","meth","total")
  
  CCC_auto <- dt[context == "CCC" & !chr%in%c("chrL", "chrM", "chrX", "chrY"), .(meth, total)]
  CCC_frac <- CCC_auto[,sum(meth)]/CCC_auto[,sum(total)]
  
  CG_auto <- dt[substr(context, start = 1, stop = 2) == "CG"  & !chr%in%c("chrL", "chrM", "chrX", "chrY"), .(meth, total)]
  CG_frac <- CG_auto[,sum(meth)]/CG_auto[,sum(total)]
  
  CH_auto <- dt[substr(context, start = 1, stop = 2) != "CG"  & !chr%in%c("chrL", "chrM", "chrX", "chrY"), .(meth, total)]
  CH_frac <- CH_auto[,sum(meth)]/CH_auto[,sum(total)]
  
  rm(dt)
  
  return(c(CCC_frac, CG_frac, CH_frac))
}

# extracted fracs from first 32 files
system.time(
  fracs_summary0 <- do.call(rbind, 
                         parallel::mclapply(paste0("data/raw_counts/raw_counts_Liu2021/GSE130553_RAW/", file_list$Name[1:100]),
                                            calFracs,
                                            mc.cores = 16
                                            )
                         )
)

fracs_summary <- as.data.table(cbind(1:nrow(fracs_summary0), round(fracs_summary0, 5)))
colnames(fracs_summary) <- c("fileNumber","CCC_Frac", "CG_Frac", "CH_Frac")

# ==== match summary data to metadata ====
## find cell numbers that corresponds to the first 32 files
merged_summary <- fuzzyjoin::difference_left_join(x = fracs_summary, y = md_1A226[,1:5], 
                                                  by = c("CCC_Frac", "CG_Frac", "CH_Frac"), 
                                                  max_dist = 0.00002
                                                  )
colnames(merged_summary)[2:4] <- c("CCC_Frac.file", "CG_Frac.file", "CH_Frac.file")
colnames(merged_summary)[7:9] <- c("CCC_Frac.cell", "CG_Frac.cell", "CH_Frac.cell")
write.csv(merged_summary, "data/metadata/GSE130553_100cell_fracs_merged.csv")



