source("code/SETPATHS.R")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(annotatr)
library(clusterProfiler)

write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/misc/"
if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/misc/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

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

# ---- Utils ----
## Download gene annotation
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

# ---- Full case study ----
read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
res_ful.gr <- list(
  'vseq' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs"))),
  'vseq_cr' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs"))),
  'scbs' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs"))),
  'smwd' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs")))
)
annot_vseq_ful <- getGeneAnnot(res_ful.gr[['vseq']])

# ---- Homogeneous subset ----
read_dir <- 'data/interim/case_studies/luo2017mice_subset_hom/result_summary/'
res_hom.gr <- list(
  'vseq' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs"))),
  'vseq_cr' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs"))),
  'scbs' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs"))),
  'smwd' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs")))
)
annot_vseq_hom <- getGeneAnnot(res_hom.gr[['vseq']])

# ---- Comparison ----

## Compare genes from full study and homogeneous pilot study
genes_ful <- unique(annot_vseq_ful$annot.gene_id)
genes_hom <- unique(annot_vseq_hom$annot.gene_id)
sum(genes_ful %in% genes_hom) / length(genes_ful) # = 0.473794
sum(genes_hom %in% genes_ful) / length(genes_hom) # = 0.8166786

## GO analysis on non-overlap genes in two studies
method <- 'vseq'
for (main_study in c('ful', 'hom')) {
  
  if (main_study == 'ful') {
    genes_nonolap <- genes_ful[! genes_ful %in% genes_hom]
    name_seg <- 'fulNoHom'
  } else if (main_study == 'hom') {
    genes_nonolap <- genes_hom[! genes_hom %in% genes_ful]
    name_seg <- 'homNoFul'
  } else {
    stop('error')
  }
  
  for (ont in c("MF", "CC", "BP")) {
    enrich_go <- enrichGO(
      gene          = genes_nonolap,
      universe      = unique(annt_genes$gene_id),
      keyType       = "ENTREZID",
      OrgDb         = org.Mm.eg.db,
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.01,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
    
    enrich_go.df <- enrich_go %>% as.data.frame() %>% dplyr::select(- geneID)
    
    if (nrow(enrich_go.df)) {
      fwrite(enrich_go.df, 
             paste0(write_dir, 'GOenrich_ORA_ont', ont, '_vseq_', name_seg, '.csv'))
      
      myXlim <- function(ont) switch (ont,
                                      'MF' = xlim(0, 0.06),
                                      'CC' = xlim(0, 0.06),
                                      'BP' = xlim(0, 0.06))
      
      dotplot(enrich_go, showCategory = 30) + 
        scale_color_gradient(low = "darkblue", high = "yellow") + 
        myXlim(ont) +
        ggtitle(paste0("Gene Ratio in GO Analysis - ", methodName(method)),
                subtitle = paste0("Ontology Category: ", ontologyName(ont))) +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))
      ggsave(paste0(plot_dir, 'GOenrich_ORA_ont', ont, '_vseq_', name_seg, '.png'),
             height = 12, width = 7)
    }
  }
}



# ---- Compare genes from my studies to marker genes in Luo 2017 paper ----
marker_genes <- read.csv("data/metadata/metadata_luo2017/luo_supp_tables/NIHMS893063-supplement-Table_S3_csv.csv", skip = 1)
idx <- grep("Cluster", colnames(marker_genes))
colnames(marker_genes)[idx] <- marker_genes[1,idx]
marker_genes <- marker_genes[-1, ]

sum(genes_ful %in% marker_genes$Gene) / length(genes_ful) # = 0.07840956
sum(genes_hom %in% marker_genes$Gene) / length(genes_hom) # = 0.07740235

## Number of unique genes mapped to each VMR
annot_vseq_ful %>% 
  select(seqnames, start, end, annot.symbol) %>%
  filter(!is.na(annot.symbol)) %>%
  group_by(seqnames, start, end) %>%
  summarise(n_genes = length(unique(annot.symbol))) %>%
  arrange(desc(n_genes))
# # A tibble: 20,055 × 4
# # Groups:   seqnames, start [20,055]
#    seqnames     start       end n_genes
#    <fct>        <int>     <int>   <int>
#  1 chr18     37832028  37832824      25
#  2 chr18     37781536  37782240      22
#  3 chr18     37787001  37787460      22
#  4 chr18     37162831  37166399      17
#  5 chr18     37180194  37182725      17
#  6 chr18     37099793  37102818      16
#  7 chr18     37026173  37027213      15
#  8 chr18     37052965  37055313      15
#  9 chr12    109711840 109712252      10
# 10 chr12    109582771 109588735       9
# # … with 20,045 more rows
# # ℹ Use `print(n = ...)` to see more rows
