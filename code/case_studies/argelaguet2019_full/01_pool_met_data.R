source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")

read_dir_met <- here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation_formatted', 'met_formatted')
write_dir <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', 'met', 'vmrseq', 'input')
if (!file.exists(write_dir)) dir.create(write_dir)

# Individual cell files
cell_file_names <- list.files(read_dir_met)

# Read in metadata
md <- fread(here('data', 'metadata', 'argelaguet2019', 'argelaguet2019_full_met&rna_sample_metadata_processed.csv'))

# Sanity check whether all cells in metadata exist
stopifnot(all(na.omit(md$file_met) %in% cell_file_names))

md$cell_file_dirs <- here(read_dir_met, md$file_met)

## ---- Process methylation data ----
# Pool data
# data.pool(cellFiles = md$cell_file_dirs, writeDir = write_dir, chrNames = 1:2, colData = md)
# data.pool(cellFiles = md$cell_file_dirs, writeDir = write_dir, chrNames = 3:5, colData = md)
# data.pool(cellFiles = md$cell_file_dirs, writeDir = write_dir, chrNames = 6:9, colData = md)
# data.pool(cellFiles = md$cell_file_dirs, writeDir = write_dir, chrNames = 10:14, colData = md)
data.pool(cellFiles = md$cell_file_dirs, writeDir = write_dir, chrNames = 15:19, colData = md)


# ---- Examine how cell cycle is related to proportion of intermediate methylation levels ----

computePropItmd <- function(cellFile) {
  df <- fread(cellFile)
  prop <- sum(df$meth_read > 0 & df$meth_read < df$total_read) / nrow(df)
  return(prop)
}

md$prop_itmd_met <- do.call(c, parallel::mclapply(md$cell_file_dirs, computePropItmd, mc.cores = 8))

md %>%
  ggplot(aes(Phase, prop_itmd_met)) +
  geom_violin() +
  ylim(0, 0.05) +
  ylab('Proportion of CpGs with intermediate methylation levels')
md %>%
  ggplot(aes(Phase=='S', prop_itmd_met)) +
  geom_violin() +
  ylim(0, 0.05) +
  ylab('Proportion of CpGs with intermediate methylation levels')

md %>%
  ggplot(aes(S.Score, prop_itmd_met, color = Phase)) +
  geom_point(size = 0.5) +
  ylim(0, 0.05) +
  ylab('Proportion of CpGs with intermediate methylation levels')

md %>%
  ggplot(aes(G2M.Score, prop_itmd_met, color = Phase)) +
  geom_point(size = 0.5) +
  ylim(0, 0.05) +
  ylab('Proportion of CpGs with intermediate methylation levels')
