source('code/SETPATHS.R')
devtools::load_all('../vmrseq-package/vmrseq/')
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(scales)

read_dir_met <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', '04_summarize_output_met')
read_dir_rna <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', 'rna')
write_dir <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', '06_met_rna_corr')
if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- here('plots', 'case_studies', 'argelaguet2019_full', '06_met_rna_corr')
if (!file.exists(plot_dir)) dir.create(plot_dir)

# ---- Fixed arguments ----
max_dist <- 1000
corr_method <- 'spearman'
upstream <- 2000; downstream <- 0
top_pct_var_genes <- 0.1 # Percentage of top variable genes

go_p_cutoff <- 0.01 # adjusted p-val cutoff in GO enrichment analysis
go_q_cutoff <- 0.05 # adjusted q-val cutoff in GO enrichment analysis

# ---- Load data ----

md <- fread(here('data', 'metadata', 'argelaguet2019', 'argelaguet2019_full_met&rna_sample_metadata_processed.csv'))

vmr.se <- loadHDF5SummarizedExperiment(here(read_dir_met, 'vmrseq_regionSummary_vmrs'))
colnames(vmr.se) <- md$sample
colData(vmr.se) <- cbind(colData(vmr.se), md)

promoter.se <- loadHDF5SummarizedExperiment(here(read_dir_met, 'promoters_up2kb_regionSummary'))
colnames(promoter.se) <- md$sample
dim(promoter.se) # 14933   939

rna_0.se <- readRDS(here(read_dir_rna, 'SingleCellExperiment_rna_vst.rds'))
dim(rna.se) #  15559   939

## Only take top variable genes
var_rna <- rowVars(assays(rna_0.se)$vst_counts)
rna.se <- rna_0.se[order(var_rna, decreasing = T)[1:(top_pct_var_genes*nrow(rna_0.se))], ]

## Subset to genes with present promoter methylation
promoter.se <- promoter.se[match(granges(rna.se)$symbol, granges(promoter.se)$symbol) %>% na.omit, ]
rna.se      <- rna.se     [match(granges(promoter.se)$symbol, granges(rna.se)$symbol) %>% na.omit, ] 
dim(promoter.se) # 1493   939
dim(rna.se)      # 1493   939
all(granges(rna.se)$symbol == granges(promoter.se)$symbol) # TRUE

## Make colnames of vmr.se, promoter.se and rna.se align
all(colnames(rna.se) == md$id_rna) # TRUE
colnames(rna.se) <- md$sample 

quantile(width(granges(vmr.se)))
# 0%    25%    50%    75%   100% 
# 13   1226   2724   6670 223003 
quantile(width(granges(promoter.se)))
#   0%  25%  50%  75% 100% 
# 2000 2000 2000 2000 2000  

sum(countOverlaps(granges(promoter.se), granges(vmr.se)) > 0) / nrow(promoter.se) 
# ~57% promoters overlap with VMRs

# # ---- GO analysis on the selected top variable gene ----
# #### (Run in singularity) ####
# enrich_go <- clusterProfiler::enrichGO(
#   gene          = unique(rowData(rna.se)$symbol),
#   universe      = unique(rowData(rna_0.se)$symbol),
#   keyType       = 'SYMBOL',
#   OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
#   ont           = 'ALL',
#   pool          = FALSE,
#   pAdjustMethod = 'BH',
#   pvalueCutoff  = go_p_cutoff,
#   qvalueCutoff  = go_q_cutoff,
#   readable      = TRUE
# )
# saveRDS(enrich_go, here(write_dir, paste0('enrichGO_ontALL_topVarGenesInRNA.rds')))
# ####
# 
# 
# ## Plot GO enrichment analysis results
# enrich_go <- readRDS(here(write_dir, paste0('enrichGO_ontALL_topVarGenesInRNA.rds')))
# 
# clusterProfiler::dotplot(enrich_go,
#                          showCategory = 50,
#                          label_format = 50) +
#   scale_color_gradient(low = 'darkblue', high = 'yellow') +
#   ggtitle(paste0('Gene Ratio in GO Analysis'),
#           subtitle = paste0('Ontology Category: BP, MF, CC')) +
#   theme(plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5))
# ggsave(here(plot_dir, paste0('enrichGO_dotplot_ontALL_topVarGenesInRNA.png')), width = 7, height = 12)

# ---- VMRs as ref: Compare VMRs and promoters in terms of correlation with gene expression ----

## Use VMRs as reference: pick the closest gene to each VMR
near_vmr.hits <- distanceToNearest(granges(vmr.se), granges(rna.se)) %>%
  as.data.frame() %>%
  filter(distance <= max_dist)
vmr_near.se      <- vmr.se     [near_vmr.hits$queryHits, ]
rna_near.se      <- rna.se     [near_vmr.hits$subjectHits, ]
promoter_near.se <- promoter.se[near_vmr.hits$subjectHits, ]

## Infer genomic context of these VMRs
vmr_context <- rep(paste0('VMRs outside gene but within ', max_dist, '-bp distance'), nrow(vmr_near.se))
vmr_context[which(granges(vmr_near.se) %over% granges(rna_near.se))] <- 'VMRs overlapping with gene body'
vmr_context[which(granges(vmr_near.se) %over% granges(promoter_near.se))] <- 'VMRs overlapping with promoter'
table(vmr_context)

## Compute correlation
vmr.mtx <- (assays(vmr_near.se)$M / assays(vmr_near.se)$Cov) %>% as.matrix
rna.mtx <- (assays(rna_near.se)$vst_counts)
promoter.mtx <- (assays(promoter_near.se)$M / assays(promoter_near.se)$Cov) %>% as.matrix

N <- nrow(vmr.mtx)
corr_vmr      <- map_dbl(.x = 1:N, 
                         .f = ~ cor(vmr.mtx[.x, ], rna.mtx[.x, ], method = corr_method, use = 'pairwise.complete.obs'))
corr_promoter <- map_dbl(.x = 1:N, 
                         .f = ~ cor(promoter.mtx[.x, ], rna.mtx[.x, ], method = corr_method, use = 'pairwise.complete.obs'))

n_avail_cell_vmr      <- map_dbl(.x = 1:N, 
                                 .f = ~ sum(!is.na(vmr.mtx[.x, ])      & !is.na(rna.mtx[.x, ])))
n_avail_cell_promoter <- map_dbl(.x = 1:N, 
                                 .f = ~ sum(!is.na(promoter.mtx[.x, ]) & !is.na(rna.mtx[.x, ])))

rowData(vmr_near.se) <- cbind(
  rowData(vmr_near.se),
  data.frame(
    rna_gene_name          = granges(rna_near.se)$symbol,
    rna_gene_start         = start(rna_near.se),
    rna_gene_end           = end(rna_near.se),
    rna_gene_strand        = strand(rna_near.se),
    vmr_to_gene_dist       = near_vmr.hits$distance,
    vmr_rna_corr           = corr_vmr, 
    vmr_n_avail_cell       = n_avail_cell_vmr,
    loglik_diff            = granges(vmr_near.se)$loglik_diff,
    promoter_rna_corr      = corr_promoter,
    promoter_n_avail_cell  = n_avail_cell_promoter,
    context                = vmr_context
  )
)

saveRDS(granges(vmr_near.se), here(write_dir, paste0('GRanges_correlation_', corr_method, '_metNrna_vmrsAsRef_maxDist', max_dist, 'bp_promotorUp2kb.rds')))
saveHDF5SummarizedExperiment(vmr_near.se, here(write_dir, paste0('SummarizedExperiment_correlation_', corr_method, '_metNrna_vmrsAsRef_maxDist', max_dist, 'bp_promotorUp2kb')),
                             replace = T)

# ---- GO analysis on high corr genes ----
vmr_corr.df <- readRDS(here(write_dir, paste0('GRanges_correlation_', corr_method, '_metNrna_vmrsAsRef_maxDist', max_dist, 'bp_promotorUp2kb.rds'))) %>% 
  as.data.frame %>%
  arrange(desc(loglik_diff))

# QC on VMRs and promoters
plot.df <- vmr_corr.df %>%
  filter(vmr_n_avail_cell >= 10 & promoter_n_avail_cell >= 10) ## ensure min 10 available cells covered

##### (!! run in singularity) ###########
# Genes with high positive/negative corr with VMRs
enrichGOonHighlightGenes <- function(object, cutoff, direction) {
  
  corr_column <- paste0(object, '_rna_corr')
  if (direction == 'positive') {
    hightlight.df <- plot.df %>% filter(plot.df[corr_column] > cutoff) 
  } else if (direction == 'negative') {
    hightlight.df <- plot.df %>% filter(plot.df[corr_column] < cutoff) 
  }
  
  hightlight.df <- hightlight.df %>% 
    arrange(desc(abs(hightlight.df[corr_column]))) %>%
    filter(!duplicated(paste0(rna_gene_name, context)))
  fwrite(hightlight.df, here(write_dir, paste0('hightlighGenes_', object, 'CorrWithRNA_', direction, '_cutoff', abs(cutoff),'_promotorUp2kb.txt')))
  
  ## GO analysis on Genes with high positive corr with VMRs but not promoters
  for (ont in c('CC', 'BP', 'MF')) {
    enrich_go <- clusterProfiler::enrichGO(
      gene          = unique(hightlight.df$rna_gene_name),
      universe      = unique(rowData(rna.se)$symbol),
      keyType       = 'SYMBOL',
      OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
      ont           = ont,
      pAdjustMethod = 'BH',
      pvalueCutoff  = go_p_cutoff,
      qvalueCutoff  = go_q_cutoff,
      readable      = TRUE
    )
    saveRDS(enrich_go, here(write_dir, paste0('enrichGO_ont', ont, '_', object, 'CorrWithRNA_', direction, '_cutoff', abs(cutoff), '_promotorUp2kb.rds')))
    
    n_signif_go <- sum(enrich_go@result$p.adjust < enrich_go@pvalueCutoff)
    if (n_signif_go > 0) {
      ## Dotplot 
      clusterProfiler::dotplot(enrich_go,
                               showCategory = n_signif_go,
                               label_format = 80) +
        scale_color_gradient(low = 'darkblue', high = 'yellow') +
        ggtitle(paste0('Gene Ratio in GO Analysis'),
                subtitle = paste0('Ontology Category: ', ont)) +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))
      ggsave(here(plot_dir, paste0('enrichGO_dotplot_ont', ont, '_', object, 'CorrWithRNA_', direction, '_cutoff', abs(cutoff), '_promotorUp2kb.png')),
             height = 9, width = 9)
    }
  }
}

enrichGOonHighlightGenes(object = 'vmr', cutoff = 0.2, direction = 'positive')
enrichGOonHighlightGenes(object = 'vmr', cutoff = 0.3, direction = 'positive')
enrichGOonHighlightGenes(object = 'vmr', cutoff = -0.2, direction = 'negative')
enrichGOonHighlightGenes(object = 'vmr', cutoff = -0.3, direction = 'negative')
enrichGOonHighlightGenes(object = 'promoter', cutoff = -0.2, direction = 'negative')
enrichGOonHighlightGenes(object = 'promoter', cutoff = -0.3, direction = 'negative')
####################


# ---- Interpretation plots on GO results ----
getGenePresenceMatrix <- function(enrich_go) {
  signif_go <- enrich_go@result %>% filter(p.adjust < go_p_cutoff) 
  gene_list <- signif_go %>% pull(geneID) %>% strsplit('/') %>% unlist() %>% unique()
  mat <- matrix(nrow = nrow(signif_go), ncol = length(gene_list),
                dimnames = list(signif_go$Description, gene_list))
  for (i in 1:nrow(mat)) {
    gene_in_go <- signif_go$geneID[i] %>% strsplit('/') %>% unlist()
    mat[i, ] <- as.integer(colnames(mat) %in% gene_in_go)
  }
  return(mat)
}


plotGOonHighlightGenes <- function(ont, object, cutoff, direction, n_clust) {
  
  hightlight.df <- fread(here(write_dir, paste0('hightlighGenes_', object, 'CorrWithRNA_', direction, '_cutoff', abs(cutoff),'_promotorUp2kb.txt')))
  
  enrich_go <- readRDS(here(write_dir, paste0('enrichGO_ont', ont, '_', object, 'CorrWithRNA_', direction, '_cutoff', abs(cutoff), '_promotorUp2kb.rds')))
  n_signif_go <- sum(enrich_go@result$p.adjust < enrich_go@pvalueCutoff)
  
  if (n_signif_go > 2) {

    ## Corresponded clustering from treeplot above
    enrich_go_2 <- enrichplot::pairwise_termsim(enrich_go, showCategory = n_signif_go)
    go_corr.mat <- enrichplot:::fill_termsim(enrich_go_2, keep = 1:n_signif_go)
    go_hc <- stats::hclust(stats::as.dist(1 - go_corr.mat), method = 'ward.D')
    # go_clust.df <- stats::cutree(go_hc, k = n_clust) %>% as.factor %>% as.data.frame
    go_clust.df <- stats::cutree(go_hc, h = 1) %>% as.factor %>% as.data.frame
    colnames(go_clust.df) <- 'GO term cluster'
    n_clust <- length(unique(go_clust.df$`GO term cluster`))
    
    annt_color <- wesanderson::wes_palette("Cavalcanti1")[1:n_clust]
    names(annt_color) <- 1:n_clust
    
    ## Matrix of gene presence indicator in each enriched GO term
    go_gene.mat <- getGenePresenceMatrix(enrich_go)
    graphics.off()
    png(here(plot_dir, paste0('enrichGO_heatmap_ont', ont, '_', object, 'CorrWithRNA_', direction, '_cutoff', abs(cutoff), '_promotorUp2kb.png')),
        width = 1200, height = 500)
    pheatmap::pheatmap(go_gene.mat,
                       color = c('darkgrey', 'white'),
                       border_color = 'white',
                       cellwidth = 13,
                       cellheight = 13,
                       legend = F,
                       cluster_rows = go_hc,
                       cutree_rows = n_clust,
                       cluster_cols = T,
                       clustering_method = 'ward.D2',
                       treeheight_col = 0,
                       show_colnames = T,
                       show_rownames = T,
                       angle_col = '45',
                       annotation_row = go_clust.df,
                       annotation_colors = list(`GO term cluster` = annt_color),
                       annotation_legend = F)
    dev.off()
  }
}

plotGOonHighlightGenes(ont = 'MF', object = 'vmr', cutoff = 0.2, direction = 'positive', n_clust = 2)
plotGOonHighlightGenes(ont = 'BP', object = 'vmr', cutoff = 0.2, direction = 'positive', n_clust = 5)
plotGOonHighlightGenes(ont = 'CC', object = 'vmr', cutoff = 0.2, direction = 'positive', n_clust = 2) # no output since there are only 2 significant GO terms
plotGOonHighlightGenes(ont = 'MF', object = 'vmr', cutoff = 0.3, direction = 'positive', n_clust = 2)
plotGOonHighlightGenes(ont = 'BP', object = 'vmr', cutoff = 0.3, direction = 'positive', n_clust = 4)



# ---- Plot correlation comparison (highlighting gene with high positive correlation with VMRs) ----
vmr_corr.df <- readRDS(here(write_dir, paste0('GRanges_correlation_', corr_method, '_metNrna_vmrsAsRef_maxDist', max_dist, 'bp_promotorUp2kb.rds'))) %>% 
  as.data.frame %>%
  arrange(desc(loglik_diff))

# QC on VMRs and promoters
plot.df <- vmr_corr.df %>%
  filter(vmr_n_avail_cell >= 10 & promoter_n_avail_cell >= 10) ## ensure min 10 available cells covered

cutoff <- 0.2
enrich_go_BP <- readRDS(here(write_dir, paste0('enrichGO_ontBP_vmrCorrWithRNA_positive_cutoff', abs(cutoff), '_promotorUp2kb.rds')))
gene_in_BP <- enrich_go_BP@result %>% filter(p.adjust < go_p_cutoff) %>% pull(geneID) %>% strsplit('/') %>% unlist() %>% unique()
enrich_go_MF <- readRDS(here(write_dir, paste0('enrichGO_ontMF_vmrCorrWithRNA_positive_cutoff', abs(cutoff), '_promotorUp2kb.rds')))
gene_in_MF <- enrich_go_MF@result %>% filter(p.adjust < go_p_cutoff) %>% pull(geneID) %>% strsplit('/') %>% unlist() %>% unique()

hightlight.df <- fread(here(write_dir, paste0('hightlighGenes_vmrCorrWithRNA_positive_cutoff', abs(cutoff), '_promotorUp2kb.txt'))) %>%
  mutate(in_go = rna_gene_name %in% c(gene_in_BP, gene_in_MF))

bg_color_1 <- 'steelblue'
bg_color_2 <- 'darkgrey'
highlight_color <- 'black'
text_annt.df <- plot.df %>%
  group_by(context) %>%
  summarise(
    pct_dots_in_bg_1 = sum(
      (vmr_rna_corr < 0) & (promoter_rna_corr < 0) & (abs(vmr_rna_corr) > abs(promoter_rna_corr))
    ) / length(vmr_rna_corr),
    pct_dots_in_bg_2 = sum(
      (vmr_rna_corr < 0) & (promoter_rna_corr < 0) & (abs(vmr_rna_corr) < abs(promoter_rna_corr))
    ) / length(vmr_rna_corr)
  ) %>%
  mutate(text_bg_1 = paste0(round(pct_dots_in_bg_1*100, 1), '% dots in blue area'),
         text_bg_2 = paste0(round(pct_dots_in_bg_2*100, 1), '% dots in grey area'))

plot.df %>%
  ggplot(aes(promoter_rna_corr, vmr_rna_corr)) +
  geom_polygon(data = data.frame(x = c(-1, 0, 0), y = c(-1, -1, 0)), 
               aes(x, y), fill = bg_color_1, alpha = 0.2) +
  geom_polygon(data = data.frame(x = c(-1, -1, 0), y = c(-1, 0, 0)), 
               aes(x, y), fill = bg_color_2, alpha = 0.2) +
  geom_rug(aes(color = context), alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = 'grey', linetype = 3) +
  geom_hline(yintercept = 0, color = 'grey', linetype = 3) +
  geom_vline(xintercept = 0, color = 'grey', linetype = 3) +
  geom_hline(yintercept = cutoff, color = highlight_color, linetype = 'dashed') +
  geom_point(aes(color = context), size = 1, alpha = 1) +
  geom_point(data = hightlight.df, color = highlight_color, size = 1.5, shape = 1) +
  ggrepel::geom_label_repel(data = hightlight.df,
                            aes(label = rna_gene_name, fill = in_go),
                            color         = highlight_color,
                            size          = 2,
                            box.padding   = 0.3,
                            point.padding = 0,
                            nudge_y       = 0.2,
                            max.overlaps  = 50,
                            segment.color = 'grey50') +
  geom_text(data = text_annt.df, aes(label = text_bg_1), x = 0.5, y = -0.7, color = bg_color_1) +
  geom_text(data = text_annt.df, aes(label = text_bg_2), x = 0.5, y = -0.8, color = bg_color_2) +
  facet_wrap(~ context) +
  scale_color_manual(values = c('#CC6666', '#9999CC', '#66CC99')) +
  scale_fill_manual(values = c('white', wesanderson::wes_palette("Moonrise3")[5])) +
  scale_fill_manual(values = c('white', 'lightgrey')) + # color for gene name boxes
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) +
  ylab('Spearman correlation of gene expression with VMR methylation') +
  xlab('Spearman correlation of gene expression with promoter methylation') +
  theme_classic() +
  theme(legend.position = 'none',
        panel.spacing.x = unit(c(0.1, 0.1), 'null'))
ggsave(here(plot_dir, paste0('correlation_', corr_method, '_metNrna_vmrsAsRef_maxDist', max_dist, 'bp_cutoff', abs(cutoff), '_promotorUp2kb.png')), width = 15, height = 5.5)



