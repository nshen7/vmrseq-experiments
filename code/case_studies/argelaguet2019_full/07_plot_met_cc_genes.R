source("code/SETPATHS.R")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(scales)
library(pdftools)

read_dir_met <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', '04_summarize_output_met')
read_dir_rna <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', 'rna')
write_dir <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', '07_plot_met_cc_genes')
if (!file.exists(write_dir)) dir.create(write_dir)

# ---- Fixed arguments ----
max_dist <- 0
min_n_cell <- 5

# ---- Load in data ----

# md <- fread(here('data', 'metadata', 'argelaguet2019', 'argelaguet2019_full_met&rna_sample_metadata_processed.csv')) %>%
#   mutate(Phase = factor(Phase, levels = c('G1', 'S', 'G2M')))
md <- fread(here('data', 'interim', 'case_studies', 'argelaguet2019_full', '05_plot_global_met', 'metadata_umap_regional_methyl_vseq_seed2010.txt.gz')) %>%
  mutate(Phase = factor(Phase, levels = c('G1', 'S', 'G2M')))

vmr.se <- loadHDF5SummarizedExperiment(here(read_dir_met, "vmrseq_regionSummary_vmrs"))
colnames(vmr.se) <- md$sample
colData(vmr.se) <- cbind(colData(vmr.se), md)

# promoter.se <- loadHDF5SummarizedExperiment(here(read_dir_met, "promoters_regionSummary"))
# colnames(promoter.se) <- md$sample
# dim(promoter.se) # 20525   939

rna.se <- readRDS(here(read_dir_rna, 'SingleCellExperiment_rna_vst.rds'))
dim(rna.se) #  15559   939

# ---- Sanity check that CC genes are more variable than housekeeping genes ----
plot_dir <- here('plots', 'case_studies', 'argelaguet2019_full', '07_plot_met_cc_genes')
if (!file.exists(plot_dir)) dir.create(plot_dir)

load(here('data', 'metadata', 'housekeeping_genes', 'Housekeeping_Genes_Mouse.RData'))
hk_genes.df <- Mouse_HK_genes; rm(Mouse_HK_genes)
cc_genes_1.df <- fread(here('data', 'metadata', 'reCAT', 'supp_table_1_processed.csv'))
cc_genes_2.df <- fread(here('data', 'metadata', 'reCAT', 'supp_table_2_processed.csv'))
# cc_genes_3 <- c(
#   fread(here('data', 'metadata', 'argelaguet2019', 'mouse_cell_cycle_genes_fromSeurat_S.csv')) %>% unlist,
#   fread(here('data', 'metadata', 'argelaguet2019', 'mouse_cell_cycle_genes_fromSeurat_G2M.csv')) %>% unlist
# )

hk_var <- assays(rna.se)$vst_count[granges(rna.se)$symbol %in% hk_genes.df$Gene, ] %>% rowVars
cc_1_var <- assays(rna.se)$vst_count[granges(rna.se)$symbol %in% cc_genes_1.df$GeneMouse, ] %>% rowVars
cc_2_var <- assays(rna.se)$vst_count[granges(rna.se)$symbol %in% cc_genes_2.df$GeneMouse, ] %>% rowVars
df <- rbind(
  data.frame(var = hk_var, gene_type = 'Housekeeping'),
  data.frame(var = cc_1_var, gene_type = 'Cell Cycle reCAT Supp Table 1'),
  data.frame(var = cc_2_var, gene_type = 'Cell Cycle reCAT Supp Table 2')
)
df %>%
  ggplot(aes(gene_type, var)) + 
  geom_jitter(color = 'pink', alpha = 0.4) +
  geom_boxplot(alpha = 0.5) +
  ylim(0, 10) +
  theme_classic()
ggsave(here(plot_dir, 'boxplot_geneVar_cc&hk.png'), width = 6, height = 4)

# # ---- (Optional) Use more strict CC phase definition ----
# 
# md <- md %>%
#   filter((Phase == 'S' & S.Score > 0.1 & G2M.Score < S.Score - 0.1) | 
#            (Phase == 'G2M' & G2M.Score > 0.1 & G2M.Score > S.Score + 0.1) | 
#            (Phase == 'G1' & S.Score < -0.1 & G2M.Score < 0.1))
# 
# vmr.se <- vmr.se[, colnames(vmr.se) %in% md$sample]
# colnames(vmr.se) <- md$sample
# colData(vmr.se) <- cbind(colData(vmr.se), md)
# 

# ---- From reCAT (Liu et al.) Supp Table 2: 15 high confidence cell cycle genes selected according to published literatures ----
plot_dir <- here('plots', 'case_studies', 'argelaguet2019_full', '07_plot_met_cc_genes', 'recat_supp_table_2')
if (!file.exists(plot_dir)) dir.create(plot_dir)

## (Running in Singularity!!) Process supp table 2 from reCAT paper
# supp_file <- pdf_text(here('data', 'metadata', 'reCAT', 'reCAT_supp_info.pdf')) %>% readr::read_lines()
# cc_genes_0.df <- supp_file[(grep('Supplementary Table 2', supp_file) + 2):(grep('Supplementary Table 3', supp_file) - 1)] %>%
#   str_squish() %>%
#   strsplit(split = ' ')
# cc_genes.df <- do.call(rbind, cc_genes_0.df[-1]) %>%
#   as.data.frame() %>%
#   mutate(V2 = paste(V2, V3)) %>%
#   select(-V3)
# colnames(cc_genes.df) <- cc_genes_0.df[[1]]
#
# cc_genes.df <- cc_genes.df %>% # Get corresponding mouse genes
#   mutate(GeneMouse = gprofiler2::gorth(Gene, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name)
# fwrite(cc_genes.df, here('data', 'metadata', 'reCAT', 'supp_table_2_processed.csv'))

## VMRs near these 15 genes
cc_genes.df <- fread(here('data', 'metadata', 'reCAT', 'supp_table_2_processed.csv')) %>%
  filter(GeneMouse %in% rowData(rna.se)$symbol)
rna_cc.se <- rna.se[match(cc_genes.df$GeneMouse, rowData(rna.se)$symbol),]
dim(rna_cc.se) #  15 939

vmr_cc.hits <- distanceToNearest(granges(vmr.se), granges(rna_cc.se)) %>%
  as.data.frame() %>%
  filter(distance <= max_dist)
vmr_cc.se <- vmr.se[vmr_cc.hits$queryHits, ]
rowData(vmr_cc.se) <- cbind(rowData(vmr_cc.se), cc_genes.df[vmr_cc.hits$subjectHits] %>% select(Peakstage, GeneMouse))

## Subset cells to those that are not globally hypomethylated
vmr_cc_sub.se <- vmr_cc.se[, which(colData(vmr_cc.se)$umap_1 < 0)]
vmr_cc_sub.mtx <- (assays(vmr_cc_sub.se)$M / assays(vmr_cc_sub.se)$Cov) %>% as.matrix

## Plot VMR individually
for (i in 1:nrow(vmr_cc_sub.se)) {
  df <- data.frame(colData(vmr_cc_sub.se)[c('Phase')], regional_methyl = vmr_cc_sub.mtx[i, ])
  vmr_info <- granges(vmr_cc_sub.se)[i]
  df %>%
    ggplot(aes(Phase, regional_methyl)) +
    geom_violin() +
    geom_jitter(alpha = 0.5, color = 'pink', height = 0) + 
    ggtitle(paste0('VMR ', i, ' (chr', seqnames(vmr_info), ':', start(vmr_info), '-', end(vmr_info), ')'), 
            subtitle = paste0('Gene: ', vmr_info$GeneMouse, '; Peakstage: ', vmr_info$Peakstage)) +
    xlab('Cell cycle phase') + 
    ylab('Regional methylation of VMRs') +
    theme_classic()
  ggsave(here(plot_dir, paste0('violin_regionalMethyl_vs_cellCyclePhase_individualVMR_VMR',i,'.png')))
}

smr <- data.frame(colData(vmr_cc_sub.se)[c('Phase')], t(vmr_cc_sub.mtx)) %>%
  pivot_longer(cols = -1, names_to = 'vmr_no', values_to = 'mf') %>%
  mutate(vmr_no = gsub('X(.*)', '\\1', vmr_no) %>% as.integer) %>%
  group_by(Phase, vmr_no) %>%
  summarise(n_cell  = sum(!is.na(mf)),
            mean_mf   = mean(mf, na.rm = T),
            q1_mf     = quantile(mf, 0.25, na.rm = T),
            q3_mf     = quantile(mf, 0.75, na.rm = T)) %>%
  mutate(gene_symbol    = rowData(vmr_cc.se)$GeneMouse[vmr_no],
         gene_peak_stage = rowData(vmr_cc.se)$Peakstage[vmr_no])

smr %>%
  filter(n_cell > min_n_cell) %>%
  ggplot(aes(Phase, mean_mf, color = gene_peak_stage)) +
  geom_boxplot()

smr %>%
  filter(n_cell > min_n_cell) %>%
  ggplot(aes(gene_peak_stage, mean_mf, color = Phase)) +
  geom_boxplot()

smr %>%
  group_by(vmr_no) %>%
  filter(min(n_cell) > 10, gene_peak_stage == 'G1') %>%
  ggplot(aes(Phase, mean_mf, color = gene_symbol)) +
  geom_path(aes(group = vmr_no))
smr %>%
  group_by(vmr_no) %>%
  filter(min(n_cell) > 10, gene_peak_stage == 'G2') %>%
  ggplot(aes(Phase, mean_mf, color = gene_symbol)) +
  geom_path(aes(group = vmr_no))
smr %>%
  group_by(vmr_no) %>%
  filter(min(n_cell) > 10, gene_peak_stage == 'S') %>%
  ggplot(aes(Phase, mean_mf, color = gene_symbol)) +
  geom_path(aes(group = vmr_no))




# ---- From reCAT (Liu et al.) Supp Table 1: The top 20 genes for each of the cell cycle stages (G1, S, G2/M) from Cyclebase. ---- 
plot_dir <- here('plots', 'case_studies', 'argelaguet2019_full', '07_plot_met_cc_genes', 'recat_supp_table_1')
if (!file.exists(plot_dir)) dir.create(plot_dir)

## (Running in Singularity!!) Process supp table 1 from reCAT paper
# supp_file <- pdf_text(here('data', 'metadata', 'reCAT', 'reCAT_supp_info.pdf')) %>% readr::read_lines()
# cc_genes_0.df <- supp_file[(grep('Supplementary Table 1', supp_file) + 2):(grep('Supplementary Table 2', supp_file) - 1)] %>%
#   str_squish() %>%
#   strsplit(split = ' ') 
# cc_genes.df <- do.call(rbind, cc_genes_0.df[-1]) %>%
#   as.data.frame() %>%
#   mutate(V2 = paste(V2, V3)) %>%
#   select(-V3)
# colnames(cc_genes.df) <- cc_genes_0.df[[1]]
# 
# # Get corresponding mouse genes
# orth <- gprofiler2::gorth(cc_genes.df$Gene, source_organism = "hsapiens", target_organism = "mmusculus") %>%
#   mutate(Gene = input, GeneMouse = ortholog_name) %>%
#   select(Gene, GeneMouse) 
# cc_genes.df <- cc_genes.df %>% 
#   left_join(orth) %>%
#   filter(!is.na(GeneMouse))
# fwrite(cc_genes.df, here('data', 'metadata', 'reCAT', 'supp_table_1_processed.csv'))

## VMRs near these CC genes
cc_genes.df <- fread(here('data', 'metadata', 'reCAT', 'supp_table_1_processed.csv')) %>%
  filter(GeneMouse %in% rowData(rna.se)$symbol) %>%
  mutate(Peakstage = factor(Peakstage, c('G1', 'S', 'G2')))
rna_cc.se <- rna.se[match(cc_genes.df$GeneMouse, rowData(rna.se)$symbol),]
dim(rna_cc.se) #  55 939

# Include all VMRs within range of 1000bp of CC genes
vmr_cc.hits <- distanceToNearest(granges(vmr.se), granges(rna_cc.se)) %>%
  as.data.frame() %>%
  filter(distance <= max_dist)
vmr_cc.se <- vmr.se[vmr_cc.hits$queryHits, ]
rowData(vmr_cc.se) <- cbind(rowData(vmr_cc.se), cc_genes.df[vmr_cc.hits$subjectHits, ] %>% select(Peakstage, GeneMouse))
quantile(rowData(vmr_cc.se)$num_cpg)

# # Pick 1 VMR for each gene
# vmr_cc.hits <- distanceToNearest(granges(rna_cc.se), granges(vmr.se)) %>%
#   as.data.frame() %>%
#   filter(distance <= max_dist)
# vmr_cc.se <- vmr.se[vmr_cc.hits$subjectHits, ]
# rowData(vmr_cc.se) <- cbind(rowData(vmr_cc.se), cc_genes.df[vmr_cc.hits$queryHits, ] %>% select(Peakstage, GeneMouse))

# Select VMRs
# vmr_cc.se <- vmr_cc.se[which(rowData(vmr_cc.se)$num_cpg > 10), ]

## Subset cells to those that are not globally hypomethylated
vmr_cc_sub.se <- vmr_cc.se[, which(colData(vmr_cc.se)$umap_1 < 0)]
vmr_cc_sub.mtx <- (assays(vmr_cc_sub.se)$M / assays(vmr_cc_sub.se)$Cov) %>% as.matrix


# ## Plot VMR individually
# which(granges(vmr_cc_sub.se)$pi < 0.6)
# # [1]   6  33  50  62  84  95 100 108
# i <- 33
# granges(vmr_cc_sub.se)[i]
# df <- data.frame(colData(vmr_cc_sub.se)[c('Phase')], regional_methyl = vmr_cc_sub.mtx[i, ])
# df %>%
#   ggplot(aes(Phase, regional_methyl)) +
#   geom_boxplot()


## Summary for each VMR
smr <- data.frame(colData(vmr_cc_sub.se)[c('Phase')], t(vmr_cc_sub.mtx)) %>%
  pivot_longer(cols = -1, names_to = 'vmr_no', values_to = 'mf') %>%
  mutate(vmr_no = gsub('X(.*)', '\\1', vmr_no) %>% as.integer) %>%
  group_by(Phase, vmr_no) %>%
  summarise(n_cell  = sum(!is.na(mf)),
            mean_mf = mean(mf, na.rm = T),
            sd_mf   = sd(mf, na.rm = T),
            q1_mf   = quantile(mf, 0.25, na.rm = T),
            q3_mf   = quantile(mf, 0.75, na.rm = T)) %>%
  mutate(gene_symbol    = rowData(vmr_cc.se)$GeneMouse[vmr_no],
         gene_peak_stage = rowData(vmr_cc.se)$Peakstage[vmr_no])


smr %>%
  group_by(vmr_no) %>%
  filter(n_cell > min_n_cell) %>%
  ggplot(aes(gene_peak_stage, mean_mf, color = Phase)) +
  geom_boxplot() +
  xlab('Gene peak stage') + 
  ylab('Regional methylation of VMRs averaged over cells') +
  scale_color_discrete(name = 'Cell cycle phase') +
  theme_classic()
ggsave(here(plot_dir, 'boxplot_regionalMethyl_vs_genePeakStage_coloredByPhase.png'), width = 5, height = 4)

for (this_stage in levels(smr$gene_peak_stage)) {
  smr %>%
    group_by(vmr_no) %>%
    filter(min(n_cell) > min_n_cell, gene_peak_stage == this_stage) %>%
    ggplot(aes(Phase, mean_mf, color = gene_symbol, group = vmr_no)) +
    geom_path() +
    # geom_errorbar(aes(ymin = mean_mf - sd_mf, ymax = mean_mf + sd_mf)) +
    ggtitle(paste0('Genes with peak stage ', this_stage)) +
    xlab('Cells grouped by cell cycle phase') + 
    ylab('Regional methylation of VMRs averaged over cells') +
    theme_classic()
  ggsave(here(plot_dir, paste0('path_regionalMethyl_vs_cellCyclePhase_', this_stage, '.png')), width = 5, height = 4)
}


# ---- From Seurat CC genes ---- 

## VMRs near these CC genes
m.s.genes <- fread(here('data', 'metadata', 'argelaguet2019', 'mouse_cell_cycle_genes_fromSeurat_S.csv'), header = F) %>% unlist
m.g2m.genes <- fread(here('data', 'metadata', 'argelaguet2019', 'mouse_cell_cycle_genes_fromSeurat_G2M.csv'), header = F) %>% unlist
cc_genes.df <- data.frame(
  GeneMouse = c(m.s.genes, m.g2m.genes),
  Peakstage = c(rep('S', length(m.s.genes)), rep('G2M', length(m.g2m.genes)))
) %>%
  filter(GeneMouse %in% rowData(rna.se)$symbol)

rna_cc.se <- rna.se[match(cc_genes.df$GeneMouse, rowData(rna.se)$symbol),]
dim(rna_cc.se) #  92 939

vmr_cc.hits <- distanceToNearest(granges(vmr.se), granges(rna_cc.se)) %>%
  as.data.frame() %>%
  filter(distance <= max_dist)
vmr_cc.se <- vmr.se[vmr_cc.hits$queryHits, ]
rowData(vmr_cc.se) <- cbind(rowData(vmr_cc.se), cc_genes.df[vmr_cc.hits$subjectHits,] %>% select(GeneMouse, Peakstage))

vmr_cc_sub.se <- vmr_cc.se[, which(colData(vmr_cc.se)$stage == 'E7.5')]
vmr_cc_sub.mtx <- (assays(vmr_cc_sub.se)$M / assays(vmr_cc_sub.se)$Cov) %>% as.matrix

smr <- data.frame(colData(vmr_cc_sub.se)[c('Phase')], t(vmr_cc_sub.mtx)) %>%
  pivot_longer(cols = -1, names_to = 'vmr_no', values_to = 'mf') %>%
  mutate(vmr_no = gsub('X(.*)', '\\1', vmr_no) %>% as.integer) %>%
  group_by(Phase, vmr_no) %>%
  summarise(n_cell  = sum(!is.na(mf)),
            mean_mf = mean(mf, na.rm = T),
            q1_mf   = quantile(mf, 0.25, na.rm = T),
            q3_mf   = quantile(mf, 0.75, na.rm = T)) %>%
  mutate(gene_symbol    = rowData(vmr_cc.se)$GeneMouse[vmr_no],
         gene_peak_stage = rowData(vmr_cc.se)$Peakstage[vmr_no])

smr %>%
  filter(n_cell > min_n_cell) %>%
  ggplot(aes(Phase, mean_mf, color = gene_peak_stage)) +
  geom_boxplot()

smr %>%
  group_by(vmr_no) %>%
  filter(min(n_cell) > min_n_cell, gene_peak_stage == 'G2M') %>%
  ggplot(aes(Phase, mean_mf, color = gene_symbol)) +
  geom_path(aes(group = vmr_no))
smr %>%
  group_by(vmr_no) %>%
  filter(min(n_cell) > min_n_cell, gene_peak_stage == 'S') %>%
  ggplot(aes(Phase, mean_mf, color = gene_symbol)) +
  geom_path(aes(group = vmr_no))






# ---- Heatmap with VMRs near CC genes ----
plot_dir <- here('plots', 'case_studies', 'argelaguet2019_full', '07_plot_met_cc_genes')
if (!file.exists(plot_dir)) dir.create(plot_dir)

cc_genes_1 <- fread(here('data', 'metadata', 'reCAT', 'supp_table_2_processed.csv')) %>%
  filter(GeneMouse %in% rowData(rna.se)$symbol) %>%
  pull(GeneMouse)
cc_genes_2 <- fread(here('data', 'metadata', 'reCAT', 'supp_table_1_processed.csv')) %>%
  filter(GeneMouse %in% rowData(rna.se)$symbol) %>%
  mutate(Peakstage = factor(Peakstage, c('G1', 'S', 'G2'))) %>%
  pull(GeneMouse)
cc_genes_3 <- rbind(fread(here('data', 'metadata', 'argelaguet2019', 'mouse_cell_cycle_genes_fromSeurat_S.csv'), header = F),
                    fread(here('data', 'metadata', 'argelaguet2019', 'mouse_cell_cycle_genes_fromSeurat_G2M.csv'), header = F)) %>%
  unlist()

cc_genes <- c(cc_genes_1, cc_genes_2, cc_genes_3)
# cc_genes <- c(cc_genes_3)
rna_cc.gr <- granges(rna.se[match(cc_genes, rowData(rna.se)$symbol) %>% na.omit(),])
length(rna_cc.gr) #  162

vmr_cc.hits <- distanceToNearest(granges(vmr.se), rna_cc.gr) %>%
  as.data.frame() %>%
  filter(distance <= max_dist)
vmr_cc.se <- vmr.se[vmr_cc.hits$queryHits, ]
rowData(vmr_cc.se) <- cbind(rowData(vmr_cc.se), cc_genes.df[vmr_cc.hits$subjectHits,] %>% select(GeneMouse, Peakstage))

plotRegionalMethylHeatmap <- function(stage = NULL, dissim_metric = "manhattan", hclust_metric = "ward", umap_metric = "euclidean") {
  
  ## Subset to one developmental stage
  if (!is.null(stage)) vmr_cc.se <- vmr_cc.se[, which(colData(vmr_cc.se)$stage == stage)]
  name_seg <- ifelse(is.null(stage), yes = '', no = paste0('_stage', stage))
  
  PHASE_COLOR <- c("#0072B2", "#ABDDA4", "#CC79A7")
  names(PHASE_COLOR) <- c('G1', 'G2M', 'S')
  
  mtx <- (assays(vmr_cc.se)$M / assays(vmr_cc.se)$Cov) %>% as.matrix
  md <- colData(vmr_cc.se) %>% as.data.frame
  
  # Sort matrix and annotation in `Phase` order
  idx <- order(md$Phase)
  md_sorted <- md[idx, ]; rownames(md_sorted) <- md_sorted$sample
  mtx_sorted <- mtx[,idx]; colnames(mtx_sorted) <- md_sorted$sample
  
  # Get dissimilarity matrix for columns (i.e., cells)
  d_mat_cols <- cluster::daisy(t(mtx_sorted), metric = dissim_metric, stand = FALSE) %>% as.matrix()
  rownames(d_mat_cols) <- colnames(d_mat_cols) <- md_sorted$sample
  
  graphics.off() 
  png(here(plot_dir, paste0("heatmap_cell_dissimilarity_ccVMRs", name_seg, ".png")), height = 750, width = 900)
  heatmap_palette <- rev(colorRampPalette(RColorBrewer::brewer.pal(8, name = "YlOrRd"))(21))
  pheatmap::pheatmap(d_mat_cols,
                     # cluster_rows = cluster_rows,
                     # cluster_cols = cluster_results$hclust_obj,
                     # treeheight_row = 0,
                     border_color = NA,
                     color = heatmap_palette,
                     cluster_cols = T,
                     cluster_rows = T,
                     show_colnames = F,
                     show_rownames = F,
                     na_col = "grey90",
                     annotation_col = md_sorted %>% select(Phase, stage, lineage10x_2),
                     annotation_colors = list(Phase = PHASE_COLOR),
                     main = paste0("Cell-to-cell dissimilarity heatmap using VMRs near CC genes")
  )
  dev.off()
}

plotRegionalMethylHeatmap()
plotRegionalMethylHeatmap(stage = 'E4.5')
plotRegionalMethylHeatmap(stage = 'E5.5')
plotRegionalMethylHeatmap(stage = 'E6.5')
plotRegionalMethylHeatmap(stage = 'E7.5')

