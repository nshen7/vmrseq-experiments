source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")

read_dir_rna <- here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation', 'rna')
write_dir_rna <- here('data', 'interim', 'case_studies', 'argelaguet2019', 'rna')
if (!file.exists(write_dir_rna)) dir.create(write_dir_rna)

write_dir_md <- here('data', 'metadata', 'argelaguet2019')
if (!file.exists(write_dir_md)) dir.create(write_dir_md)

# Metadata
md <- fread(here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation', 'sample_metadata.txt')) %>%
  mutate(file_met = ifelse(is.na(id_met), NA, paste0(id_met, '.tsv.gz'))) %>%
  filter(pass_metQC & pass_rnaQC & !is.na(lineage10x_2))

fwrite(md, here(write_dir_md, 'Argelaguet2019_met&rna_sample_metadata_processed.csv'))

## ---- Process gene data ----

rna_count.se <- readRDS(here(read_dir_rna, 'SingleCellExperiment.rds'))
dim(rna_count.se) # 18345  2482
colnames(rna_count.se)

# Subset to only the samples with QC-passed met&rna profiles 
rna_count_sub.se <- rna_count.se[, md$id_rna]
dim(rna_count_sub.se) # 18345   939

saveRDS(rna_count_sub.se, here(write_dir_rna, 'SingleCellExperiment_met&rna.rds'))
