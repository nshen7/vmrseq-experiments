source("code/SETPATHS.R")
library(doParallel)
registerDoParallel(8)  

read_dir <- here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation')
write_dir <- here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation_formatted')
if(!file.exists(write_dir)) dir.create(write_dir)

## ---- Load in data ----
read_dir_met <- here(read_dir, 'met', 'cpg_level')

# File names of cells with methylome profiled
met_files <- list.files(read_dir_met)

# Metadata
md <- fread(here(read_dir, 'sample_metadata.txt')) %>%
  mutate(file_met = ifelse(is.na(id_met), NA, paste0(id_met, '.tsv.gz')))

# Checking if all cells are available
all(na.omit(md$file_met) %in% met_files) # = TRUE

## ---- Process CpG files into suitable input file for vmrseq::data.pool ----
write_dir_met <- here(write_dir, 'met_formatted')
if(!file.exists(write_dir_met)) dir.create(write_dir_met)

formatCellFiles <- function(cell_file_name) {
  
  # Read in CpG files of a cell
  cell.df <- fread(here(read_dir_met, cell_file_name))
  if (!all(colnames(cell.df) == c('chr', 'pos', 'met_reads', 'nonmet_reads', 'rate')))
    stop('Colnames does not align with expectations!')
  
  # Format into vmrseq::data.pool input 
  formatted.df <- with(cell.df, data.frame(chr        = chr,
                                           pos        = pos,
                                           strand     = '*',
                                           meth_read  = met_reads,
                                           total_read = met_reads + nonmet_reads))
  
  # Write out formatted CpG files and gzip it
  write_dir_file <- here(write_dir_met, gsub('(.*).gz', '\\1', cell_file_name))
  write_tsv(formatted.df, write_dir_file)
  R.utils::gzip(write_dir_file, remove = T, overwrite = T)
  
}

# Run (took several mins)
foreach (file_name = met_files) %dopar% formatCellFiles(file_name)


