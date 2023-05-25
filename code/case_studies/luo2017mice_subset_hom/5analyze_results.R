source('code/SETPATHS.R')
setwd('/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/')
devtools::load_all('../vmrseq-package/vmrseq/')
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(scales)
library(pheatmap)
library(annotatr)
library(clusterProfiler)

read_dir <- 'data/interim/case_studies/luo2017mice_subset_hom/result_summary/'
plot_dir <- 'plots/case_studies/luo2017mice_subset_hom/comparison/'
if (!file.exists(plot_dir)) dir.create(plot_dir)

md <- fread('data/interim/case_studies/luo2017mice_subset_hom/metadata_luo2017mice_subset_hom.csv')
rownames(md) <- md$sample

# ==== read in results from various methods ====
sites.gr <- readRDS(paste0(read_dir, 'cpg_sites.rds'))
res <- list(
  vseq = loadHDF5SummarizedExperiment(paste0(read_dir, 'vmrseq_regionSummary_vmrs')),
  vseq_cr = loadHDF5SummarizedExperiment(paste0(read_dir, 'vmrseq_regionSummary_crs')),
  scbs = loadHDF5SummarizedExperiment(paste0(read_dir, 'scbs_regionSummary_vmrs')),
  smwd = loadHDF5SummarizedExperiment(paste0(read_dir, 'smallwood_regionSummary_vmrs'))
)
res.gr <- map(res, granges)

methods <- c('vmrseq', 'CR in vmrseq', 'scbs', 'smallwood')
methods <- factor(methods, levels = methods)

################################
###### General exploration #####
################################

# Color settings
COLORS <- RColorBrewer::brewer.pal(n = 6, name = 'PuOr')[-3]
COLORVALUES <- c('vmrseq' = COLORS[1], 'CR in vmrseq' = COLORS[2], 
                 'scbs' = COLORS[3], 'smallwood' = COLORS[4])

# Number & percentage of sites in detected regions
n_sites <- sapply(res_region, function(se) findOverlaps(GenomicRanges::reduce(granges(se)), sites.gr) %>% length())
ggplot() + 
  geom_bar(aes(x = methods, y = n_sites, color = methods, fill = methods), stat = 'identity') + 
  scale_y_continuous(labels = scientific, limits = c(0, 2.5e6), name = 'Total number of CpGs in detected regions') +
  scale_color_manual(values = COLORVALUES) + 
  scale_fill_manual(values = COLORVALUES) + 
  xlab('Methods') + 
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid = element_line(color = 'light grey'))
ggsave(paste0(plot_dir, 'barplot_nSites_vs_methods.png'), height = 4, width = 4)  

# Percentage of sites in detected regions
pct <- n_sites / length(sites.gr)
ggplot() + 
  geom_bar(aes(x = methods, y = pct, color = methods, fill = methods), stat = 'identity') + 
  scale_y_continuous(labels = percent, limits = c(0, 0.118), name = 'Percentage of CpGs in detected regions') +
  scale_color_manual(values = COLORVALUES) + 
  scale_fill_manual(values = COLORVALUES) + 
  xlab('Methods') + 
  theme(legend.position = 'none',
        panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid = element_line(color = 'light grey'))
ggsave(paste0(plot_dir, 'barplot_pctSites_vs_methods.png'), height = 4, width = 4)  


# Distribution of CpG number per detected region
n_sites_avail <- lapply(res_region, function(se) countOverlaps(granges(se), sites.gr))
ggplot() +
  geom_density(aes(n_sites_avail[[1]], color = names(COLORVALUES)[1])) + 
  geom_density(aes(n_sites_avail[[2]], color = names(COLORVALUES)[2])) + 
  geom_density(aes(n_sites_avail[[3]], color = names(COLORVALUES)[3])) + 
  geom_density(aes(n_sites_avail[[4]], color = names(COLORVALUES)[4])) + 
  scale_x_continuous(limits = c(0, 400), name = 'Number of CpG per detected region') +
  theme(panel.background = element_rect(fill = 'white', color = 'black'), panel.grid = element_line(color = 'light grey')) +
  scale_color_manual(name = 'Methods', breaks = names(COLORVALUES), values = COLORVALUES)
ggsave(paste0(plot_dir, 'densityplot_nCpGsPerRegion.png'), height = 4, width = 7)  


# Number of covered regions relative to total number of detected regions per cell
n_regn_avail <- lapply(res_region, function(se) colSums(assays(se)$Cov > 0)/nrow(se))
ggplot() +
  geom_density(aes(n_regn_avail[[1]], color = names(COLORVALUES)[1])) + 
  geom_density(aes(n_regn_avail[[2]], color = names(COLORVALUES)[2])) + 
  geom_density(aes(n_regn_avail[[3]], color = names(COLORVALUES)[3])) + 
  geom_density(aes(n_regn_avail[[4]], color = names(COLORVALUES)[4])) + 
  scale_x_continuous(labels = percent, limits = c(0, 1), name = 'Percentage of covered regions per cell') +
  theme(panel.background = element_rect(fill = 'white', color = 'black'), panel.grid = element_line(color = 'light grey')) +
  scale_color_manual(name = 'Methods', breaks = names(COLORVALUES), values = COLORVALUES)
ggsave(paste0(plot_dir, 'densityplot_pctCoveredRegionPerCell.png'), height = 4, width = 7)  


###############################################################
###### Heatmap of methylation in highly variable regions ######
###############################################################

heatmapTopRegions <- function(n_top, method, dissim_metric = 'manhattan', hclust_metric = 'ward') {
  # `n_top`: number of top regions
  # `method`: either 'vseq', 'scbs', 'smwd' or 'scmet'  ('vseq' not available since no ranking available)
  
  # Extact top VMRs 
  metric <- switch(method,
                   'vseq' = granges(res_region[[method]])$loglik_diff,
                   'scbs' = granges(res_region[[method]])$mcols.var,
                   'smwd' = granges(res_region[[method]])$var_lb,
                   'scmet' = granges(res_region[[method]])$tail_prob)
  method_name <- switch(method,
                        'vseq' = 'vmrseq',
                        'scbs' = 'scbs',
                        'smwd' = 'Smallwood',
                        'scmet' = 'scMET')
  
  top_ind <- order(metric, decreasing = TRUE)[1:n_top]
  top.se <- res_region[[method]][top_ind]
  top_MF <- as.matrix(assays(top.se)$M / assays(top.se)$Cov)
  
  # Sort matrix and annotation in `Neuron_type3` order
  idx <- order(md$Neuron_type3)
  md_sorted <- md[idx, ]; rownames(md_sorted) <- md_sorted$sample
  top_MF_sorted <- top_MF[,idx]; colnames(top_MF_sorted) <- md_sorted$sample
  
  # Get hierarchical clustering on rows and columns
  d_mat_rows <- as.matrix(cluster::daisy(top_MF_sorted, metric = dissim_metric, stand = FALSE))
  d_mat_rows[is.na(d_mat_rows)] <- mean(d_mat_rows, na.rm = TRUE)
  cluster_rows <- cluster::agnes(d_mat_rows, diss = TRUE, method = hclust_metric)
  d_mat_cols <- as.matrix(cluster::daisy(t(top_MF_sorted), metric = dissim_metric, stand = FALSE))
  d_mat_cols[is.na(d_mat_cols)] <- mean(d_mat_cols, na.rm = TRUE)
  cluster_cols <- cluster::agnes(d_mat_cols, diss = TRUE, method = hclust_metric)
  
  graphics.off() 
  png(paste0(plot_dir, 'heatmap_top', n_top, 'regions_', method, '.png'), height = 900, width = 900)
  heatmap_palette <- colorRampPalette(RColorBrewer::brewer.pal(8, name = 'YlOrRd'))(21)
  pheatmap(top_MF_sorted[cluster_rows$order, cluster_cols$order],
           # cluster_rows = cluster_rows,
           # cluster_cols = cluster_results$hclust_obj,
           # treeheight_row = 0,
           border_color = NA,
           color = heatmap_palette,
           cluster_cols = F,
           cluster_rows = F,
           show_colnames = F,
           show_rownames = F,
           na_col = 'grey90',
           annotation_col = md_sorted[cluster_cols$order,] %>% 
             mutate(names = rownames(md_sorted)[cluster_cols$order]) %>%
             column_to_rownames('names') %>% 
             select(Neuron_type3, Neuron_type1),
           main = paste0('Regional Mean Methylation of Top ', n_top, ' Regions from ', method_name)
  )
  dev.off()
} 

heatmapTopRegions(n_top = 100, method = 'vseq')
heatmapTopRegions(n_top = 100, method = 'scbs')
heatmapTopRegions(n_top = 100, method = 'smwd') 

heatmapTopRegions(n_top = 300, method = 'vseq')
heatmapTopRegions(n_top = 300, method = 'scbs')
heatmapTopRegions(n_top = 300, method = 'smwd') 

heatmapTopRegions(n_top = 1000, method = 'vseq')
heatmapTopRegions(n_top = 1000, method = 'scbs')
heatmapTopRegions(n_top = 1000, method = 'smwd') 









######################################
###### Gene enrichment analysis ######
######################################
write_dir <- 'data/interim/case_studies/luo2017mice_subset_hom/gene_enrichment/'
if (!file.exists(write_dir)) dir.create(write_dir)

## Download gene annotation
types_gene <- grep("mm10_genes.*", builtin_annotations(), value = T)
types_gene <- types_gene[types_gene != "mm10_genes_intergenic"]
annt_genes <- build_annotations(genome = 'mm10', annotations = types_gene)

## Util functions
methodName <- function(method) switch (method,
                                       'vseq' = 'vmrseq',
                                       'vseq_cr' = 'vmrseq CRs',
                                       'scbs' = 'scbs',
                                       'smwd' = 'Smallwood',
                                       'scmet' = 'scMET')
ontologyName <- function(ont) switch(ont,
                                     'MF' = 'Molecular Function',
                                     'CC' = 'Cellular Component',
                                     'BP' = 'Biological Process')

getGeneAnnot <- function(regions.gr){
  method_annt <- annotate_regions(
    regions = regions.gr,
    annotations = annt_genes,
    ignore.strand = TRUE,
    quiet = FALSE
  ) %>% data.frame()
  return(method_annt)
}

performGO <- function(method) {
  regions.gr <- res.gr[[method]]
  method_annt <- getGeneAnnot(regions.gr)
  for (ont in c('MF', 'CC', 'BP')) {
    enrich_go <- enrichGO(
      gene          = unique(method_annt$annot.gene_id),
      universe      = unique(annt_genes$gene_id),
      keyType       = 'ENTREZID',
      OrgDb         = org.Mm.eg.db,
      ont           = ont,
      pAdjustMethod = 'BH',
      pvalueCutoff  = 0.01,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
    
    fwrite(enrich_go %>% as.data.frame() %>% dplyr::select(- geneID), 
           paste0(write_dir, 'GOenrich_ORA_ont', ont, '_', method, '.csv'))
    
    myXlim <- function(ont) switch (ont,
                                    'MF' = xlim(0, 0.05),
                                    'CC' = xlim(0, 0.05),
                                    'BP' = xlim(0, 0.05))
    
    dotplot(enrich_go, showCategory = 30) + 
      scale_color_gradient(low = 'darkblue', high = 'yellow') + 
      myXlim(ont) +
      ggtitle(paste0('Gene Ratio in GO Analysis - ', methodName(method)),
              subtitle = paste0('Ontology Category: ', ontologyName(ont))) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    ggsave(paste0(plot_dir, 'GOenrich_ORA_ont', ont, '_', method, '.png'),
           height = 12, width = 7)
  }
}

performGO(method = 'vseq')
performGO(method = 'vseq_cr')
performGO(method = 'scbs')
performGO(method = 'smwd')



