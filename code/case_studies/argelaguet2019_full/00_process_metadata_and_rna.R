source("code/SETPATHS.R")
library(SummarizedExperiment)

read_dir_rna <- here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation', 'rna')
write_dir_rna <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', 'rna')
if (!file.exists(write_dir_rna)) dir.create(write_dir_rna)

write_dir_md <- here('data', 'metadata', 'argelaguet2019')
if (!file.exists(write_dir_md)) dir.create(write_dir_md)

## ---- Metadata ----

md <- fread(here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation', 'sample_metadata.txt')) %>%
  mutate(file_met = ifelse(is.na(id_met), NA, paste0(id_met, '.tsv.gz'))) %>%
  filter(pass_metQC & pass_rnaQC & !is.na(lineage10x_2))

fwrite(md, here(write_dir_md, 'argelaguet2019_full_met&rna_sample_metadata_processed_0.csv'))

# ---- RNA-seq: Subset to only the samples with QC-passed met&rna profiles ----
rna_count.se <- readRDS(here(read_dir_rna, 'SingleCellExperiment.rds'))
rna_count_sub.se <- rna_count.se[, md$id_rna]
dim(rna_count_sub.se) # 18345   939

# Remove genes that are covered by less than 5 cells
rna_count_sub.se <- rna_count_sub.se[-which(rowMaxs(assays(rna_count_sub.se)$counts) < 5), ]
dim(rna_count_sub.se) # 17568   939

seqlevelsStyle(rna_count_sub.se) <- 'NCBI'

saveRDS(rna_count_sub.se, here(write_dir_rna, 'SingleCellExperiment_rna.rds'))


# ---- Tag cell with cell cycles inferred from gene expression ----
rna_count_sub.se <- readRDS(here(write_dir_rna, 'SingleCellExperiment_rna.rds'))
md <- fread(here(write_dir_md, 'argelaguet2019_full_met&rna_sample_metadata_processed_0.csv'))

## !! run in local terminal (had issues with installing Seurat in singularity) 
########
## Ref vignette: https://satijalab.org/seurat/articles/cell_cycle_vignette.html
library(Seurat)
library(gprofiler2)
m.s.genes <- gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
m.g2m.genes <- gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
fwrite(m.s.genes %>% as.data.frame, here(write_dir_md, 'mouse_cell_cycle_genes_fromSeurat_S.csv'), col.names = F)
fwrite(m.g2m.genes %>% as.data.frame, here(write_dir_md, 'mouse_cell_cycle_genes_fromSeurat_G2M.csv'), col.names = F)

# Create our Seurat object and complete the initalization steps
counts <- assays(rna_count_sub.se)$counts
rownames(counts) <- rowData(rna_count_sub.se)$symbol
counts.srt <- CreateSeuratObject(counts = counts)
counts.srt <- NormalizeData(counts.srt)
counts.srt <- CellCycleScoring(counts.srt, 
                               s.features = m.s.genes, 
                               g2m.features = m.g2m.genes, 
                               set.ident = TRUE)

all(colnames(counts.srt) == md$id_rna)
md$S.Score = counts.srt$S.Score
md$G2M.Score = counts.srt$G2M.Score
md$Phase = counts.srt$Phase
fwrite(md, here(write_dir_md, 'argelaguet2019_full_met&rna_sample_metadata_processed.csv'))
########


# ---- normalize rna-seq counts by vst (!!need to be run locally due to environment conflicts) ----

rna.se <- readRDS(here(write_dir_rna, 'SingleCellExperiment_rna.rds'))

out_vst <- sctransform::vst(assays(rna.se)$counts)
rna.se <- rna.se[match(rownames(out_vst$y), rownames(rna.se)),]
assays(rna.se)$vst_counts <- out_vst$y
saveRDS(rna.se, here(write_dir_rna, 'SingleCellExperiment_rna_vst.rds'))

