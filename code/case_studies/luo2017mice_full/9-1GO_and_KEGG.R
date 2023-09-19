source("code/SETPATHS.R")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(annotatr)
library(clusterProfiler)

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/gene_set_enrichment/"
if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/gene_set_enrichment/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

res.gr <- list(
  'vseq' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs"))),
  'vseq_cr' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs"))),
  'scbs' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs"))),
  'smwd' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "Smallwood_regionSummary_vmrs"))),
  'scmet' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "scmet_regionSummary_vmrs")))
)


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

getRegions <- function(method, top_n_regions = NULL, reduce_smwd = TRUE){
  
  regions.gr <- res.gr[[method]]
  if (method != "vseq_cr") {
    metric <- switch(method,
                     "vseq" = regions.gr$loglik_diff,
                     "scbs" = regions.gr$mcols.var,
                     "smwd" = regions.gr$var_lb,
                     "scmet" = regions.gr$gamma)
    regions.gr <- regions.gr[order(metric, decreasing = TRUE)]
  }
  
  if (!is.null(top_n_regions)) {
    if (method %in% c("vseq_cr")) stop("Method does not provide rank!")
    regions.gr <- regions.gr[1:top_n_regions]
  }
  
  if (method == "smwd" & reduce_smwd) {
    regions.gr <- GenomicRanges::reduce(regions.gr)
  }
  
  return(regions.gr)
}

# ---- How many VMRs overlaps with genes? ----

computeProportion <- function(method) {
  
  helper <- function(n, method) {
    
    if (n != Inf) {
      if (method == 'vseq_cr') 
        return(NULL)
      if (n > length(getRegions(method, NULL))) 
        return(NULL)
      regions.gr <- getRegions(method, top_n_regions = n, reduce_smwd = FALSE)
    } else {
      regions.gr <- getRegions(method, top_n_regions = NULL, reduce_smwd = FALSE) 
    }
      
    distance <- distanceToNearest(regions.gr, annt_genes) %>% data.frame() %>% pull(distance)
    pct_in_genes <- sum(distance == 0) / length(regions.gr)
    
    return(data.frame(method, top_n_regions = length(regions.gr), pct_in_genes))
  }
  
  n_all <- c(100, 1000, 3000, 10000, 30000, Inf)
  summary.df <- map_dfr(.x = n_all, .f = helper, method = method)
  return(summary.df)
} 

pct_in_genes.df <- rbind(
  computeProportion('vseq'),
  computeProportion('vseq_cr'),
  computeProportion('scbs'),
  computeProportion('smwd'),
  computeProportion('scmet')
) %>%
  mutate(method = factor(method)) %>%
  mutate(method = fct_recode(method, 
                             "vmrseq" = "vseq", "vmrseq CRs" = "vseq_cr", 
                             "scbs" = "scbs", "Smallwood" = "smwd",
                             "scMET" = "scmet")) 

# Color settings
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = COLORS[1], "vmrseq CRs" = COLORS[2],
                 "scbs" = COLORS[3], "Smallwood" = COLORS[4], "scMET" = COLORS[5])


pct_in_genes.df %>%
  ggplot(aes(top_n_regions, pct_in_genes, color = method)) + 
  geom_point() +
  geom_path() + 
  geom_point(data = pct_in_genes.df %>% filter(method == 'vmrseq CRs'), 
             aes(top_n_regions, pct_in_genes, color = method), 
             size = 3) +
  scale_color_manual(values = COLORVALUES) +
  scale_x_log10() + ylim(0, 1) +
  # scale_fill_manual(values = COLORVALUES) +
  xlab("N Top Regions") + ylab("Percentage of detected regions overlaping genes") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_line(color = "light grey"))
ggsave(here(plot_dir, '..', paste0("point_pctInGene_vs_nTopRgions.png")), height = 4, width = 5)



# ---- GO Over Representation Analysis ----
performGO <- function(method, top_n_regions = NULL) {
  
  regions.gr <- getRegions(method, top_n_regions) 
  method_annt <- getGeneAnnot(regions.gr)

  for (ont in c("MF", "CC", "BP")) {
    enrich_go <- enrichGO(
      # gene          = unique(method_annt$annot.symbol),
      # universe      = unique(annt_genes$symbol),
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
    
    seg <- ifelse(is.null(top_n_regions), yes = "", no = paste0("_top", top_n_regions, "regions"))
    fwrite(enrich_go %>% as.data.frame() %>% dplyr::select(- geneID), 
           paste0(write_dir, 'GOenrich_ORA_ont', ont, '_', method, seg, ".csv"))
    
    myXlim <- function(ont) switch (ont,
                                    'MF' = xlim(0.002, 0.035),
                                    'CC' = xlim(0.003, 0.045),
                                    'BP' = xlim(0.008, 0.045))
    
    dotplot(enrich_go, showCategory = 30) + 
      scale_color_gradient(low = "darkblue", high = "yellow") + 
      myXlim(ont) +
      ggtitle(paste0("Gene Ratio in GO Analysis - ", methodName(method)),
              subtitle = paste0("Ontology Category: ", ontologyName(ont))) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    ggsave(paste0(plot_dir, 'GOenrich_ORA_ont', ont, '_', method, seg, ".png"),
           height = 12, width = 7)
  }
}

performGO(method = "vseq")
performGO(method = "vseq_cr")
performGO(method = "scbs")
performGO(method = "smwd")
performGO(method = "scmet")

# ---- KEGG Over Representation Analysis ----

## Need to be run outside of singularity Rstudio!
performKEGG <- function(method, top_n_regions = NULL) {
  
  regions.gr <- getRegions(method, top_n_regions) # ranked regions
  method_annt <- getGeneAnnot(regions.gr)
  
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
         paste0(write_dir, 'KEGGenrich_ORA_', method, seg, ".csv"))
  
  dotplot(enrich_kegg, showCategory = 30) + 
    scale_color_gradient(low = "darkblue", high = "yellow") + 
    xlim(0.005, 0.07) +
    ggtitle(paste0("Gene Ratio in KEGG Analysis - ", methodName(method))) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(plot_dir, 'KEGGenrich_ORA_', method, seg, ".png"),
         height = 12, width = 7)
}

performKEGG(method = "vseq")
performKEGG(method = "vseq_cr")
performKEGG(method = "scbs")
performKEGG(method = "smwd")

