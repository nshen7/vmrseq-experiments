source("code/SETPATHS.R")
devtools::load_all('../vmrseq-package/vmrseq/')
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(scales)
library(pheatmap)
library(annotatr)
library(clusterProfiler)

read_dir_0 <- 'data/interim/case_studies/luo2017mice_subset_hom/result_summary/'
read_dir <- 'data/interim/case_studies/luo2017mice_subset_hom_v2/result_summary/'
plot_dir <- 'plots/case_studies/luo2017mice_subset_hom_v2/comparison/'
if (!file.exists(plot_dir)) dir.create(plot_dir)

md <- fread('data/interim/case_studies/luo2017mice_subset_hom_v2/metadata_luo2017mice_subset_hom.csv')
rownames(md) <- md$sample

# ==== read in results from various methods ====
sites.gr <- readRDS(paste0(read_dir_0, 'cpg_sites.rds'))
res <- list(
  vseq = loadHDF5SummarizedExperiment(paste0(read_dir, 'vmrseq_regionSummary_vmrs'))#,
  # vseq_cr = loadHDF5SummarizedExperiment(paste0(read_dir, 'vmrseq_regionSummary_crs')),
  # scbs = loadHDF5SummarizedExperiment(paste0(read_dir, 'scbs_regionSummary_vmrs')),
  # smwd = loadHDF5SummarizedExperiment(paste0(read_dir, 'smallwood_regionSummary_vmrs'))
)
res.gr <- map(res, granges)

methods <- c('vmrseq', 'CR in vmrseq', 'scbs', 'smallwood')
methods <- factor(methods, levels = methods)




######################################
###### Gene enrichment analysis ######
######################################
write_dir <- 'data/interim/case_studies/luo2017mice_subset_hom_v2//gene_enrichment/'
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


