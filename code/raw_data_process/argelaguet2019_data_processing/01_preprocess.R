source("code/SETPATHS.R")
library(doParallel)
registerDoParallel(8)  

read_dir <- here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation')
write_dir <- here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation_formatted')
if(!file.exists(write_dir)) dir.create(write_dir)

CHR_NUMS <- as.character(c(1:19, 'X', 'Y'))
CHR_NAMES <- paste0('chr', CHR_NUMS)

## ---- Util data and function ----
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10

# CpG coordinates in reference genome
CPG_COORDS <- map(
  CHR_NAMES, 
  ~ matchPattern('CG', genome[[.x]]) %>% as.data.frame %>% pull(start)
)
names(CPG_COORDS) <- CHR_NAMES

flipStrandAndCollapse <- function(cell.df) {
  
  if (!all(colnames(cell.df) == c('chr', 'pos', 'met_reads', 'nonmet_reads', 'rate')))
    stop('Colnames does not align with expectations!')
  
  for (i in 1:length(CHR_NUMS)) {
    # Coords of current chromosome
    ref_coords  <- CPG_COORDS[[CHR_NAMES[i]]]
    
    # Sanity check if all CpGs align with + or - strand on reference genome
    cell_coords <- cell.df %>% filter(chr == CHR_NUMS[i]) %>% pull(pos)
    if (!all(cell_coords %in% ref_coords | cell_coords %in% (ref_coords+1)))
      stop('Not all data rows are in reference genome!')
    
    # Flip strand for CpGs on the opposite strand of the ref genome
    cell.df <- cell.df %>%
      mutate(pos = ifelse(test = chr == CHR_NUMS[i] & pos %in% (ref_coords+1),
                          yes  = pos - 1,
                          no   = pos))
    # cat(i, ' ')
  }
  
  # Collapse reads on same CpG
  cell_collapsed.df <- cell.df %>%
    group_by(chr, pos) %>%
    summarise(met_reads = sum(met_reads), 
              nonmet_reads = sum(nonmet_reads))
  return(cell_collapsed.df)
}

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
  cell.df <- cell.df %>% filter(chr %in% CHR_NUMS)
  
  # Flip strand and collapse reads on the same CpGs
  cell_collapsed.df <- flipStrandAndCollapse(cell.df)
  
  # Format into vmrseq::data.pool input 
  formatted.df <- with(cell_collapsed.df, data.frame(chr        = chr,
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

# ----- misc -----
# Distribution of proportion of sites with intermediate methylation level
countItmdProportion <- function(file) {
  cell.df <- fread(here(write_dir_met, file))
  mf <- cell.df$meth_read/cell.df$total_read 
  return(sum(mf > 0 & mf < 1)/length(mf))
}
itmd_prop_argelaguet <- unlist(parallel::mclapply(list.files(write_dir_met), countItmdProportion, mc.cores = 8))
quantile(itmd_prop_argelaguet)
#          0%         25%         50%         75%        100% 
# 0.000000000 0.002869492 0.004985169 0.008636923 0.322671732 


# Distribution of proportion of sites with <=2 reads
countUnder2readProportion <- function(file) {
  cell.df <- fread(here(write_dir_met, file))
  return(sum(cell.df$total <= 2)/nrow(cell.df))
}
under2read_prop_argelaguet <- unlist(parallel::mclapply(list.files(write_dir_met), countUnder2readProportion, mc.cores = 8))
quantile(under2read_prop_argelaguet)

