# Purpose of this script:
# To investigate primitive endoderm cells because they have low global methylation levels (~=0.25)

source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
register(MulticoreParam(8))

read_dir_met <- here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation_formatted', 'met_formatted')
write_dir <- here('data', 'interim', 'case_studies', 'argelaguet2019', 'side_analysis')
if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- here('plots', 'case_studies', 'argelaguet2019', 'side_analysis')
if (!file.exists(plot_dir)) dir.create(plot_dir)

## ---- Load metadata ----
md <- fread(here('data', 'metadata', 'argelaguet2019', 'Argelaguet2019_met&rna_sample_metadata_processed.csv')) %>%
  filter(lineage10x_2 == 'Primitive_endoderm')  

all(md$stage == 'E4.5') # = TRUE

## ---- Gather methylation profiles of primitive endoderm cells ----

wrapper <- function(file) {
  df <- fread(file) %>%
    filter(chr %in% as.character(1:19)) %>%
    mutate(M = meth_read / total_read) %>% 
    filter(M == 0 | M == 1) %>% 
    select(chr, pos, M) %>%
    mutate(M = as.integer(M))
  return(df)
}

# File names of primitive endoderm cells
files <- here(read_dir_met, md$file_met[md$lineage10x_2 == 'Primitive_endoderm'])
pe.list <- lapply(files, wrapper)

## ---- Estimate transition probability from primitive endoderm cells ----

tp <- vmrseq::tp.estimate(
  list = pe.list,
  max_dist_bp = 2000,
  buffer_bp = 3000,
  lags = 1:10,
  BPPARAM = bpparam(),
  degree = 2,
  span = 0.02
)
saveRDS(tp, here(write_dir, 'transitProbs_primitive_endoderm_cell.rds'))

# plot loess-fitted transition probs
vmrseq::tp.plot(tp)
ggsave(here(plot_dir, 'transitProbs_primitive_endoderm_cell.png'), width = 8, height = 6)
