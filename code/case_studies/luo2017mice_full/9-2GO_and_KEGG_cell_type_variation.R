source("code/SETPATHS.R")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(annotatr)
library(clusterProfiler)
library(here)

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/gene_set_enrichment/"
if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/gene_set_enrichment/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

# ---- load gene annotations ----

types_gene <- grep("mm10_genes.*", builtin_annotations(), value = T)
types_gene <- types_gene[types_gene != "mm10_genes_intergenic"]
annt_genes <- build_annotations(genome = 'mm10', annotations = types_gene)

getGeneAnnot <- function(regions.gr){
  
  method_annt <- annotate_regions(
    regions = regions.gr,
    annotations = annt_genes,
    ignore.strand = TRUE,
    quiet = FALSE
  ) %>% data.frame()
  
  return(method_annt)
}

# ---- utils ----
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

## Subset VMRs to certain subset with large/small variance across cell types
subsetRegions <- function(method, cutoff, lower_tail, top_n_regions = NULL) {
  
  SE <- readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_', method, '.rds')))
  
  if (!is.null(top_n_regions)) {
    metric <- switch(method,
                     "vseq" = rowData(SE)$loglik_diff,
                     "scbs" = rowData(SE)$mcols.var,
                     "smwd" = rowData(SE)$var_lb,
                     "scmet" = rowData(SE)$gamma)
    idx <- order(metric, decreasing = T)[1:top_n_regions]
    SE <- SE[idx, ]
  }
  
  # Compute variance of regional methylation across cell types
  rowData(SE)$cell_type_var <- rowVars(as.matrix(assays(SE)$reginal_methyl), na.rm = TRUE)

  if (lower_tail) {
    regions.gr <- granges(SE) %>% subset(cell_type_var < cutoff)
  } else {
    regions.gr <- granges(SE) %>% subset(cell_type_var > cutoff)
  }
  
  return(regions.gr)
}

# ---- GO ORA on VMRs with low variance across cell types ---- 

performGO <- function(method, cutoff, lower_tail, top_n_regions = NULL) {
  
  regions.gr <- subsetRegions(method, cutoff, lower_tail, top_n_regions)
  method_annt <- getGeneAnnot(regions.gr)
  message('Found ', length(unique(method_annt$annot.gene_id)), ' annotated genes...')
  
  for (ont in c("MF", "CC", "BP")) {
    enrich_go <- enrichGO(
      gene          = unique(method_annt$annot.gene_id),
      universe      = unique(annt_genes$gene_id),
      keyType       = "ENTREZID",
      OrgDb         = org.Mm.eg.db,
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.01,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
    # enrich_go %>% as.data.frame() %>% dplyr::select(-geneID) %>% View
    
    seg <- ifelse(is.null(top_n_regions), 
                  yes = '',
                  no = paste0('_selectedFromTop', top_n_regions, 'Regions'))
    fwrite(enrich_go %>% as.data.frame() %>% dplyr::select(-geneID), 
           here(write_dir, paste0('GOenrich_ORA_ont', ont, 
                                  seg, 
                                  '_', ifelse(lower_tail, 'low', 'high'), 'Var',
                                  '_', 'cutoff', cutoff,
                                  '_', method, ".csv"))
    )
    
    myXlim <- function(ont) switch (ont,
                                    'MF' = xlim(0, 0.07),
                                    'CC' = xlim(0, 0.07),
                                    'BP' = xlim(0, 0.07))
    
    dotplot(enrich_go, showCategory = 30) + 
      scale_color_gradient(low = "darkblue", high = "yellow") + 
      myXlim(ont) +
      ggtitle(paste0("Gene Ratio in GO Analysis - ", methodName(method)),
              subtitle = paste0("Ontology Category: ", ontologyName(ont))) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    ggsave(paste0(plot_dir, 'GOenrich_ORA_ont', ont, seg, 
                  '_', ifelse(lower_tail, 'low', 'high'), 'Var',
                  '_', 'cutoff', cutoff,
                  '_', method, '.png'),
           height = 12, width = 7)
  }
}

# Picked 0.02 as cutoff in reference to 95% quantile of variance distribution in smallwoods
performGO(method = 'vseq', cutoff = 0.02, lower_tail = TRUE)
performGO(method = 'vseq_cr', cutoff = 0.02, lower_tail = TRUE)
performGO(method = 'scbs', cutoff = 0.02, lower_tail = TRUE)
performGO(method = 'smwd', cutoff = 0.02, lower_tail = TRUE)
performGO(method = 'scmet', cutoff = 0.02, lower_tail = TRUE)


# ---- (failed due to package issue) KEGG ORA on VMRs with low variance across cell types ----

## Need to be run outside of singularity Rstudio!
performKEGG <- function(method, cutoff, lower_tail, top_n_regions = NULL) {
  
  regions.gr <- subsetRegions(method, cutoff, lower_tail, top_n_regions)
  method_annt <- getGeneAnnot(regions.gr)
  message('Found ', length(unique(method_annt$annot.gene_id)), ' annotated genes...')
  
  enrich_kegg <- enrichKEGG(
    gene          = unique(method_annt$annot.gene_id),
    universe      = unique(annt_genes$gene_id),
    organism      = "mmu",
    keyType       = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    use_internal_data = FALSE
  ) 
  
  seg <- ifelse(is.null(top_n_regions), yes = "", no = paste0("_top", top_n_regions, "regions"))
  fwrite(enrich_kegg %>% as.data.frame() %>% dplyr::select(- geneID), 
         paste0(write_dir, 'KEGGenrich_ORA_', seg, 
                '_', ifelse(lower_tail, 'low', 'high'), 'Var',
                '_', 'cutoff', cutoff,
                '_', method, ".csv"))
  
  dotplot(enrich_kegg, showCategory = 30) + 
    scale_color_gradient(low = "darkblue", high = "yellow") + 
    xlim(0.005, 0.07) +
    ggtitle(paste0("Gene Ratio in KEGG Analysis - ", methodName(method))) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(plot_dir, 'KEGGenrich_ORA_', seg, 
                '_', ifelse(lower_tail, 'low', 'high'), 'Var',
                '_', 'cutoff', cutoff,
                '_', method, '.png'),
         height = 12, width = 7)
}

performKEGG(method = 'vseq', cutoff = 0.02, lower_tail = TRUE)
performKEGG(method = 'vseq_cr', cutoff = 0.02, lower_tail = TRUE)
performKEGG(method = 'scbs', cutoff = 0.02, lower_tail = TRUE)
performKEGG(method = 'smwd', cutoff = 0.02, lower_tail = TRUE)
