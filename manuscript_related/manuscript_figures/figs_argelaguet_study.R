source("code/SETPATHS.R")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
source('manuscript_related/manuscript_figures/utils.R')

read_dir <- "data/interim/case_studies/argelaguet2019_full/"
plot_dir <- "manuscript_related/manuscript_figures/argelaguet_study"
write_dir <- here("manuscript_related", "manuscript_tables")
if (!file.exists(plot_dir)) dir.create(plot_dir)

df_umap <- fread(here(read_dir, "05_plot_global_met", "metadata_umap_regional_methyl_vseq_seed2010.txt.gz"))

### Color settings
STAGECOLORS <- RColorBrewer::brewer.pal(9, name = 'YlGn')[c(1,4,6,9)]
names(STAGECOLORS) <- levels(factor(df_umap$stage))
PHASECOLORS <- RColorBrewer::brewer.pal(6, name = 'BuPu')[c(1,3,6)]
names(PHASECOLORS) <- levels(factor(df_umap$Phase))

ct <- levels(factor(df_umap$lineage10x_2))
CTCOLORS <- colorRampPalette(RColorBrewer::brewer.pal(8, name = 'Accent'))(length(ct))
names(CTCOLORS) <- ct

# ---- UMAPs ----

## For plotting purpose, adjust x coords
df_umap$umap_1[df_umap$umap_1 > 5] <- df_umap$umap_1[df_umap$umap_1 > 5] - 15

## colored by time points
df_umap %>% 
  slice_sample(n = nrow(df_umap), replace = FALSE) %>%
  ggplot(aes(umap_1, umap_2)) +
  # geom_point(aes(color = stage), size = 1.5, shape = 16, alpha = 1) +
  # scale_color_manual(values = STAGECOLORS, name = 'Developmental\nStage') +
  geom_point(aes(fill = stage), size = 1.5, shape = 21, alpha = 0.8, stroke = 0.3, color = 'grey16') +
  scale_fill_manual(values = STAGECOLORS, name = 'Developmental\nStage') +
  theme_void() +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = c(.85, .75),
        legend.title.align = 0.5)
ggsave(here(plot_dir, "umap_regional_methyl_coloredByStage.png"), width = 4, height = 4)


for (this_stage in unique(df_umap$stage)) {
  
  # colored by stage
  df_umap %>% 
    filter(stage == this_stage) %>%
    ggplot(aes(umap_1, umap_2)) +
    # geom_point(aes(color = stage), size = 1.5, shape = 16, alpha = 1) +
    # scale_color_manual(values = STAGECOLORS[this_stage], name = 'Developmental\nStage') +
    geom_point(aes(fill = stage), size = 1.5, shape = 21, alpha = 0.8, stroke = 0.3, color = 'grey16') +
    scale_fill_manual(values = STAGECOLORS[this_stage], name = 'Developmental\nStage') +
    xlim(min(df_umap$umap_1), max(df_umap$umap_1)) +
    ylim(min(df_umap$umap_2), max(df_umap$umap_2)) +
    theme_void() +
    theme(strip.text.x = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = c(.85, .75),
          legend.title.align = 0.5)
  ggsave(here(plot_dir, paste0("umap_regional_methyl_coloredByStage_stage", this_stage, ".png")), width = 4, height = 4)

  # colored by cell cycle phase
  df_umap %>% 
    filter(stage == this_stage) %>%
    ggplot(aes(umap_1, umap_2)) +
    # geom_point(aes(color = Phase), size = 1.5, shape = 16, alpha = 1) +
    # scale_color_manual(values = PHASECOLORS, name = 'Cell Cycle\nPhase') +
    geom_point(aes(fill = Phase), size = 1.5, shape = 21, alpha = 0.8, stroke = 0.3, color = 'grey16') +
    scale_fill_manual(values = PHASECOLORS, name = 'Cell Cycle\nPhase') +
    xlim(min(df_umap$umap_1), max(df_umap$umap_1)) +
    ylim(min(df_umap$umap_2), max(df_umap$umap_2)) +
    theme_void() +
    theme(strip.text.x = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = c(.85, .75),
          legend.title.align = 0.5)
  ggsave(here(plot_dir, paste0("umap_regional_methyl_coloredByPhase_stage", this_stage, ".png")), width = 4, height = 4)
  
  # colored by linearge
  df_umap %>% 
    filter(stage == this_stage) %>%
    ggplot(aes(umap_1, umap_2)) +
    # geom_point(aes(color = lineage10x_2), size = 1.5, shape = 16, alpha = 1) +
    # scale_color_manual(values = CTCOLORS, name = 'Lineage', limits = force) +
    geom_point(aes(fill = lineage10x_2), size = 1.5, shape = 21, alpha = 0.8, stroke = 0.3, color = 'grey16') +
    scale_fill_manual(values = CTCOLORS, name = 'Lineage', limits = force) +
    xlim(min(df_umap$umap_1), max(df_umap$umap_1)) +
    ylim(min(df_umap$umap_2), max(df_umap$umap_2)) +
    theme_void() +
    theme(strip.text.x = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = c(.8, .65),
          legend.title.align = 0.5)
  ggsave(here(plot_dir, paste0("umap_regional_methyl_coloredByLineage_stage", this_stage, ".png")), width = 4, height = 4)
  
}

## colored by cell cycle phase
df_umap %>% 
  ggplot(aes(umap_1, umap_2)) +
  # geom_point(aes(color = Phase), size = 1.5, shape = 16, alpha = 1) +
  # scale_color_manual(values = PHASECOLORS, name = 'Cell Cycle\nPhase') +
  geom_point(aes(fill = Phase), size = 1.5, shape = 21, alpha = 0.8, stroke = 0.3, color = 'grey16') +
  scale_fill_manual(values = PHASECOLORS, name = 'Cell Cycle\nPhase') +
  theme_void() +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = c(.85,.75),
        legend.title.align = 0.5)
ggsave(here(plot_dir, "umap_regional_methyl_coloredByPhase.png"), width = 4, height = 4)

for (this_phase in unique(df_umap$Phase)) {
  
  df_umap %>% 
    filter(Phase == this_phase) %>%
    ggplot(aes(umap_1, umap_2)) +
    geom_point(aes(fill = Phase), size = 1.5, shape = 21, alpha = 0.8, stroke = 0.3, color = 'grey16') +
    scale_fill_manual(values = PHASECOLORS[this_phase], name = 'Cell Cycle\nPhase') +
    # geom_point(aes(color = Phase), size = 1.5, shape = 16, alpha = 1) +
    # scale_color_manual(values = PHASECOLORS[this_phase], name = 'Cell Cycle\nPhase') +
    xlim(min(df_umap$umap_1), max(df_umap$umap_1)) +
    ylim(min(df_umap$umap_2), max(df_umap$umap_2)) +
    theme_void() +
    theme(strip.text.x = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          legend.position = c(.85, .75),
          legend.title.align = 0.5)
  ggsave(here(plot_dir, paste0("umap_regional_methyl_coloredByPhase_phase", this_phase, ".png")), width = 4, height = 4)
  
}

# Cell cycle phase definition
df_umap %>% 
  ggplot(aes(S.Score, G2M.Score)) + 
  geom_point(aes(fill = Phase), size = 1.5, shape = 21, color = 'grey16') +
  scale_fill_manual(values = PHASECOLORS, name = 'Cell Cycle\nPhase') +
  theme_classic() 
ggsave(here(plot_dir, paste0("point_cellCycleScore_G2M_vs_S.png")),  width = 4.5, height = 4)


# colored by linearge
df_umap %>% 
  ggplot(aes(umap_1, umap_2)) +
  geom_point(aes(fill = lineage10x_2), size = 1.5, shape = 21, alpha = 0.8, stroke = 0.3, color = 'grey16') +
  scale_fill_manual(values = CTCOLORS, name = 'Lineage', limits = force) +
  # geom_point(aes(color = lineage10x_2), size = 1.5, shape = 16, alpha = 1) +
  # scale_color_manual(values = CTCOLORS, name = 'Lineage', limits = force) +
  xlim(min(df_umap$umap_1), max(df_umap$umap_1)) +
  ylim(min(df_umap$umap_2), max(df_umap$umap_2)) +
  theme_void() +
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = c(.8, .65),
        legend.title.align = 0.5)
ggsave(here(plot_dir, paste0("umap_regional_methyl_coloredByLineage.png")), width = 4, height = 4)


# ---- Methylation and expression of example VMR-gene pairs ----

read_dir_rna <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', 'rna')
read_dir_met <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', '06_met_rna_corr')
md <- fread(here('data', 'metadata', 'argelaguet2019', 'argelaguet2019_full_met&rna_sample_metadata_processed.csv'))

rna.se <- readRDS(here(read_dir_rna, 'SingleCellExperiment_rna_vst.rds'))
colnames(rna.se) <- md$sample 

vmr_near_gene.se <- loadHDF5SummarizedExperiment(here(read_dir_met, paste0('SummarizedExperiment_correlation_spearman_metNrna_vmrsAsRef_maxDist1000bp')))
colnames(vmr_near_gene.se) <- md$sample

granges(vmr_near_gene.se) %>% 
  as.data.frame %>%
  arrange(desc(vmr_rna_corr)) %>%
  head(n = 10)
granges(vmr_near_gene.se) %>% 
  as.data.frame %>%
  arrange(vmr_rna_corr) %>%
  head(n = 10)

plotVmrGenePair <- function(positive_corr, top_i) {
  
  idx <- order(granges(vmr_near_gene.se)$vmr_rna_corr, decreasing = positive_corr)[top_i]
  message('VMR No. ', idx, ' in vmr_near_gene.se')
  this_vmr <- vmr_near_gene.se[idx,]
  this_gene_name <- granges(this_vmr)$rna_gene_name
  message('Gene symbol: ', this_gene_name)
  this_gene <- rna.se[which(granges(rna.se)$symbol == this_gene_name), ]
  message('There are ', sum(t(assays(this_vmr)$Cov > 0)), ' cells available for both methylation and RNA')
  
  df <-  data.frame(colData(this_vmr) %>% as.data.frame,
                    MF = t(assays(this_vmr)$M / assays(this_vmr)$Cov),
                    GE = t(assays(this_gene)$counts) %>% unname) %>%
    filter(!is.na(MF)) %>%
    mutate(logGE = log1p(GE))
  
  df %>%
    ggplot(aes(stage, MF)) +
    geom_jitter(height = 0.02, size = 0.3) +
    geom_violin(aes(fill = stage), alpha = 0.8, width = 1.1) +
    scale_fill_manual(values = STAGECOLORS) + 
    labs(x = "Developmental stage",
         y = "Regional average methylation") +
    theme_classic() + 
    guides(fill = 'none')
  ggsave(here(plot_dir, paste0('example_', ifelse(positive_corr, 'positive', 'negative'), 'Corr_', this_gene_name, '_MF.png')), width = 3, height = 3)
  df %>%
    ggplot(aes(stage, logGE)) +
    geom_jitter(height = 0.02, size = 0.3) +
    geom_violin(aes(fill = stage), alpha = 0.8, width = 1.1) +
    scale_fill_manual(values = STAGECOLORS) + 
    labs(x = "Developmental stage",
         y = "Log-transformed gene expression\n") +
    theme_classic() + 
    theme(axis.ticks = element_blank()) + 
    guides(fill = 'none')
  ggsave(here(plot_dir, paste0('example_', ifelse(positive_corr, 'positive', 'negative'), 'Corr_', this_gene_name, '_GE.png')), width = 3, height = 3)
}

## Picking 2 highest-ranked VMRs from top 10 correlation, on each direction
plotVmrGenePair(positive_corr = TRUE, top_i = 4) # VMR No. 3841 in vmr_near_gene.se; Gene symbol: Ccnd2
plotVmrGenePair(positive_corr = TRUE, top_i = 9) # VMR No. 5241 in vmr_near_gene.se; Gene symbol: Prtg
plotVmrGenePair(positive_corr = FALSE, top_i = 3) # VMR No. 5305 in vmr_near_gene.se; Gene symbol: Dppa5a
plotVmrGenePair(positive_corr = FALSE, top_i = 7) # VMR No. 7477 in vmr_near_gene.se; Gene symbol: Slc7a7

# ---- Correlation with genes ----
max_dist <- 1000
corr_method <- 'spearman'
upstream <- 2000; downstream <- 0
top_pct_var_genes <- 0.1 # Percentage of top variable genes

go_p_cutoff <- 0.01 # adjusted p-val cutoff in GO enrichment analysis
go_q_cutoff <- 0.05 # adjusted q-val cutoff in GO enrichment analysis

read_dir_met <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', '06_met_rna_corr')
vmr_near_gene.se <- loadHDF5SummarizedExperiment(
  here(read_dir_met, paste0('SummarizedExperiment_correlation_spearman_metNrna_vmrsAsRef_maxDist1000bp_promotorUp2kb'))
)

plot.df <- granges(vmr_near_gene.se) %>%
  as.data.frame() %>%
  arrange(desc(loglik_diff)) %>%
  mutate(context = factor(context, levels = c("VMRs overlapping with promoter", "VMRs overlapping with gene body", "VMRs outside gene but within 1000-bp distance"))) %>%
  filter(vmr_n_avail_cell >= 10 & promoter_n_avail_cell >= 10) 
fwrite(plot.df %>% select(-loglik_diff.1), here(write_dir, 'supp_table_argelaguet_gene_vmr_corr.csv'))


with(plot.df %>% filter(context == "VMRs overlapping with promoter"),
     ks.test(x = vmr_rna_corr, y = promoter_rna_corr)) # D = 0.12463, p-value = 0.01066
with(plot.df %>% filter(context == "VMRs overlapping with gene body"),
     ks.test(x = vmr_rna_corr, y = promoter_rna_corr)) # D = 0.39422, p-value < 2.2e-16
with(plot.df %>% filter(context == "VMRs outside gene but within 1000-bp distance"),
     ks.test(x = vmr_rna_corr, y = promoter_rna_corr)) # D = 0.35385, p-value = 4.984e-11

# proportion of VMRs exhibiting positive corr
plot.df %>%
  group_by(context) %>%
  summarise(prop_positive = sum(vmr_rna_corr > 0)/n())
#   context                                       prop_positive
#   <fct>                                                 <dbl>
# 1 VMRs overlapping with promoter                       0.0356
# 2 VMRs overlapping with gene body                      0.111 
# 3 VMRs outside gene but within 1000-bp distance        0.159 



# plot correlation
cutoff <- 0.2
enrich_go_BP <- readRDS(here(read_dir, '06_met_rna_corr', paste0('enrichGO_ontBP_vmrCorrWithRNA_positive_cutoff', abs(cutoff), '_promotorUp2kb.rds')))
gene_in_BP <- enrich_go_BP@result %>% filter(p.adjust < go_p_cutoff) %>% pull(geneID) %>% strsplit('/') %>% unlist() %>% unique()
enrich_go_MF <- readRDS(here(read_dir, '06_met_rna_corr', paste0('enrichGO_ontMF_vmrCorrWithRNA_positive_cutoff', abs(cutoff), '_promotorUp2kb.rds')))
gene_in_MF <- enrich_go_MF@result %>% filter(p.adjust < go_p_cutoff) %>% pull(geneID) %>% strsplit('/') %>% unlist() %>% unique()

hightlight.df <- fread(here(read_dir, '06_met_rna_corr', paste0('hightlighGenes_vmrCorrWithRNA_positive_cutoff', abs(cutoff), '_promotorUp2kb.txt'))) %>%
  mutate(context = factor(context, levels = c("VMRs overlapping with promoter", "VMRs overlapping with gene body", "VMRs outside gene but within 1000-bp distance"))) %>%
  mutate(in_go = rna_gene_name %in% c(gene_in_BP, gene_in_MF))

examples.df <- data.frame(granges(vmr_near_gene.se)) %>% 
  filter(rna_gene_name %in% c('Prtg', 'Ccnd2', 'Slc7a7', 'Dppa5a')) %>%
  mutate(context = factor(context, levels = c("VMRs overlapping with promoter", "VMRs overlapping with gene body", "VMRs outside gene but within 1000-bp distance"))) %>%
  group_by(rna_gene_name) %>%
  slice(which.max(abs(vmr_rna_corr)))

bg_color_1 <- 'steelblue'
bg_color_2 <- 'darkgrey'
highlight_color <- 'black'
text_annt.df <- plot.df %>%
  group_by(context) %>%
  summarise(
    n = n(),
    pct_dots_in_bg_1 = sum(
      (vmr_rna_corr < 0) & (promoter_rna_corr < 0) & (abs(vmr_rna_corr) > abs(promoter_rna_corr))
    ) / length(vmr_rna_corr),
    pct_dots_in_bg_2 = sum(
      (vmr_rna_corr < 0) & (promoter_rna_corr < 0) & (abs(vmr_rna_corr) < abs(promoter_rna_corr))
    ) / length(vmr_rna_corr),
    pct_dots_above_line = sum(
      (vmr_rna_corr > cutoff)
    ) / length(vmr_rna_corr),
    pct_vmr_pos_corr = sum(
      (vmr_rna_corr > cutoff)
    ) / length(vmr_rna_corr)
  ) %>% 
  mutate(text_total_n = paste0('# (VMR-gene pairs) = ', n),
         text_bg_1 = paste0(round(pct_dots_in_bg_1*100, 1), '% in blue area'),
         text_bg_2 = paste0(round(pct_dots_in_bg_2*100, 1), '% in grey area'),
         text_above_line = paste0(round(pct_dots_above_line*100, 1), '% above dashed line')) 
text_annt.df %>% select(context, pct_vmr_pos_corr)
# # A tibble: 3 × 2
#   context                                        pct_vmr_pos_corr
#   <fct>                                                     <dbl>
# 1 VMRs overlapping with promoter                           0.0119
# 2 VMRs overlapping with gene body                          0.0301
# 3 VMRs outside gene but within 1000-bp distance            0.0410

## Proportion of promoters exceeding correlation > 0.2
plot.df %>%
  filter(!duplicated(rna_gene_name)) %>%
  summarise(
    pct_promoter_pos_corr = sum(
      (promoter_rna_corr > cutoff)
    ) / length(promoter_rna_corr)
  )
#   pct_promoter_pos_corr
# 1           0.007081039


set.seed(2020)
plot.df %>%
  ggplot(aes(promoter_rna_corr, vmr_rna_corr)) +
  facet_grid(~ context) +
  geom_polygon(data = data.frame(x = c(-1, 0, 0), y = c(-1, -1, 0)), 
               aes(x, y), fill = bg_color_1, alpha = 0.2) +
  geom_polygon(data = data.frame(x = c(-1, -1, 0), y = c(-1, 0, 0)), 
               aes(x, y), fill = bg_color_2, alpha = 0.2) +
  geom_path(data = data.frame(x = c(0, 0), y = c(0, 1)), aes(x, y), color = 'grey', linetype = 3) +
  geom_path(data = data.frame(x = c(0, 1), y = c(0, 0)), aes(x, y), color = 'grey', linetype = 3) +
  geom_hline(yintercept = cutoff, color = highlight_color, linetype = 'dashed') +
  geom_point(aes(color = context), shape = 3, size = 1.5, stroke = 1, alpha = 0.6) +
  geom_point(data = examples.df, color = 'palegreen4', size = 3) +
  geom_point(data = hightlight.df, color = highlight_color, size = 3, shape = 1) +
  ggrepel::geom_label_repel(data = hightlight.df,
                            aes(label = rna_gene_name, fill = in_go),
                            color         = highlight_color,
                            size          = 2,
                            box.padding   = 0.3,
                            point.padding = 0,
                            nudge_y       = 0.2,
                            max.overlaps  = 80,
                            segment.color = 'grey50') +
  geom_text(data = text_annt.df, aes(label = text_total_n), x = 0.5, y = -0.6) +
  geom_text(data = text_annt.df, aes(label = text_bg_1), x = 0.5, y = -0.7, color = bg_color_1) +
  geom_text(data = text_annt.df, aes(label = text_bg_2), x = 0.5, y = -0.8, color = bg_color_2) +
  # geom_text(data = text_annt.df, aes(label = text_above_line), x = 0.5, y = -0.9, color = 'black') +
  geom_flat_violin_x(aes(promoter_rna_corr, -1, color = context), alpha = 0.5, width = 0.3) +
  geom_flat_violin_y(aes(-1, vmr_rna_corr,      color = context), alpha = 0.5, width = 0.3) + 
  # scale_color_manual(values = c('#CC6666', '#9999CC', '#66CC99')) +
  scale_color_manual(values = wesanderson::wes_palette("IsleofDogs1")[1:3]) +
  scale_fill_manual(values = c('white', 'lightgrey')) + # color for gene name boxes
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) +
  ylab('Spearman correlation of gene expression with VMR methylation') +
  xlab('Spearman correlation of gene expression with promoter methylation') +
  theme_classic() +
  theme(legend.position = 'none',
        panel.spacing.x = unit(c(0.05, 0.05), 'null'))
ggsave(here(plot_dir, paste0('correlation_', corr_method, '_metNrna_vmrsAsRef_maxDist', max_dist, 'bp_cutoff', abs(cutoff), '.png')), width = 18, height = 6)



# ## Same plot but without gene annotation
# plot.df %>%
#   ggplot(aes(promoter_rna_corr, vmr_rna_corr)) +
#   geom_polygon(data = data.frame(x = c(-1, 0, 0), y = c(-1, -1, 0)), 
#                aes(x, y), fill = bg_color_1, alpha = 0.2) +
#   geom_polygon(data = data.frame(x = c(-1, -1, 0), y = c(-1, 0, 0)), 
#                aes(x, y), fill = bg_color_2, alpha = 0.2) +
#   geom_path(data = data.frame(x = c(0, 0), y = c(0, 1)), aes(x, y), color = 'grey', linetype = 3) +
#   geom_path(data = data.frame(x = c(0, 1), y = c(0, 0)), aes(x, y), color = 'grey', linetype = 3) +
#   geom_point(aes(color = context), shape = 3, size = 1.5, stroke = 1, alpha = 0.6) +
#   geom_point(data = examples.df, color = 'palegreen4', size = 3) +
#   geom_text(data = text_annt.df, aes(label = text_bg_1), x = 0.5, y = -0.7, color = bg_color_1) +
#   geom_text(data = text_annt.df, aes(label = text_bg_2), x = 0.5, y = -0.8, color = bg_color_2) +
#   geom_flat_violin_x(aes(promoter_rna_corr, -1, color = context), alpha = 0.5, width = 0.3) +
#   geom_flat_violin_y(aes(-1, vmr_rna_corr,      color = context), alpha = 0.5, width = 0.3) + 
#   facet_wrap(~ context, nrow = 1) +
#   scale_color_manual(values = wesanderson::wes_palette("IsleofDogs1")[c(1,2,3)]) +
#   scale_fill_manual(values = c('white', 'lightgrey')) + # color for gene name boxes
#   scale_x_continuous(expand = c(0, 0), limits = c(-1, 1)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) +
#   ylab('Spearman correlation of gene expression with VMR methylation') +
#   xlab('Spearman correlation of gene expression with promoter methylation') +
#   theme_classic() +
#   theme(legend.position = 'none',
#         panel.spacing.x = unit(c(0.1, 0.1), 'null'))
# ggsave(here(plot_dir, paste0('correlation_', corr_method, '_metNrna_vmrsAsRef_maxDist', max_dist, 'bp_cutoff', abs(cutoff), '_noGeneAnnot.png')), width = 15, height = 5)
# 
# 
# 
# 
# # ---- global methylation level ----
# 
# read_dir_raw <- here('data', 'raw_counts', 'argelaguet2019', 'scnmt_gastrulation_formatted', 'met_formatted')
# cell_file_names <- list.files(read_dir_met)
# 
# # Read in metadata
# md <- fread(here('data', 'metadata', 'argelaguet2019', 'argelaguet2019_full_met&rna_sample_metadata_processed.csv'))
# md$cell_file_dirs <- here(read_dir_met, md$file_met)
# 
# computeGlobalMet <- function(dir) {
#   df <- fread(dir)
#   mf <- df$meth_read / df$total_read
#   return(mean(mf))
# }
# 
# md$global_met <- do.call(c, parallel::mclapply(md$cell_file_dirs, computeGlobalMet, mc.cores = 8))
# 
# # colored by linearge
# md %>% 
#   ggplot(aes(lineage10x_2, global_met)) +
#   geom_jitter(height = 0.02, size = 0.05) +
#   geom_violin(aes(fill = lineage10x_2), alpha = 0.5, width = 1) +
#   scale_fill_manual(values = CTCOLORS) + 
#   scale_y_continuous(breaks = c(0,1)) +
#   facet_grid(~ stage, scales = "free_x", space = "free_x") + 
#   guides(fill = "none") +
#   xlab('Lineage commitment') +
#   ylab('Global methylation level') +
#   theme_classic() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5))
# ggsave(here(plot_dir, paste0("violin_globalMet_vs_lineage.png")), width = 6, height = 4)
# 
