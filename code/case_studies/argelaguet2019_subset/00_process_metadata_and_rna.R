## Important:
## In this script, we subset to cell types with regular global meth level (~75%)

source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")

read_dir_rna <- here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation', 'rna')
write_dir_rna <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', 'rna')
if (!file.exists(write_dir_rna)) dir.create(write_dir_rna)

write_dir_md <- here('data', 'metadata', 'argelaguet2019')
if (!file.exists(write_dir_md)) dir.create(write_dir_md)

# Metadata
md <- fread(here(write_dir_md, 'argelaguet2019_full_met&rna_sample_metadata_processed.csv')) %>%
  filter(stage != 'E4.5' & lineage10x != 'Visceral_endoderm' & lineage10x_2 != 'Visceral_endoderm' & lineage10x_2 != 'ExE_ectoderm') ## Subsetting
fwrite(md, here(write_dir_md, 'argelaguet2019_subset_met&rna_sample_metadata_processed.csv'))

## ---- Process gene data ----

rna_count.se <- readRDS(here(read_dir_rna, 'SingleCellExperiment.rds'))
dim(rna_count.se) # 18345  2482

# Subset to only the samples with QC-passed met&rna profiles 
rna_count_sub.se <- rna_count.se[, md$id_rna]
dim(rna_count_sub.se) # 18345   726

# Remove genes that are not covered by any cell
rna_count_sub.se <- rna_count_sub.se[-which(rowMaxs(assays(rna_count_sub.se)$counts)==0), ]
dim(rna_count_sub.se) # 17265   726

seqlevelsStyle(rna_count_sub.se) <- 'NCBI'

saveRDS(rna_count_sub.se, here(write_dir_rna, 'SingleCellExperiment_rna.rds'))
