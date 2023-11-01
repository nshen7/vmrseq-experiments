source("code/SETPATHS.R")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
plot_dir <- "manuscript_related/manuscript_figures/luo_study"
if (!file.exists(plot_dir)) dir.create(plot_dir)

md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")
rownames(md) <- md$sample

res_region <- list(
  vseq = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs")),
  vseq_cr = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs")),
  scbs = loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs")),
  smwd = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs")),
  scmet = loadHDF5SummarizedExperiment(paste0(read_dir, "scmet_regionSummary_vmrs"))
)
sites.gr <- readRDS(paste0(read_dir, 'cpg_sites.rds'))

# ---- fixed arguments ----
methodName <- function(method) switch (method,
                                       'vseq' = 'vmrseq',
                                       'vseq_cr' = 'vmrseq CRs',
                                       'scbs' = 'scbs',
                                       'smwd' = 'Smallwood',
                                       'scmet' = 'scMET')
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = COLORS[1], "vmrseq CRs" = COLORS[2],
                 "scbs" = COLORS[3], "Smallwood" = COLORS[4], "scMET" = COLORS[5])

methods <- c('vmrseq', 'vmrseq CRs', 'scbs', 'Smallwood', 'scMET')
methods <- factor(methods, levels = methods)

## Colors for cell subtype
ct <- md %>%
  arrange(Neuron_type1, Neuron.type) %>%
  pull(Neuron.type) %>%
  unique()
CTCOLORS <- colorRampPalette(RColorBrewer::brewer.pal(8, name = 'Accent'))(length(ct))
names(CTCOLORS) <- ct

# Color for broad cell type
BTCOLORS <- c('darkolivegreen3', 'darkgoldenrod3')
names(BTCOLORS) <- c('Excitatory', 'Inhibitory')


# ---- Get legend for cell types ----
df_umap <- fread(here(read_dir, paste0("metadata_umap_regional_methyl_vseq_seed2010.txt.gz")))

# For UMAPs
(p <- df_umap %>% 
    ggplot() +
    # geom_point(aes(fill = Neuron.type, shape = Neuron_type1), size = 1.5, alpha = 0.8, stroke = 0.3, color = 'grey16') +
    # scale_fill_manual(values = CTCOLORS, name = "") +
    geom_point(aes(umap_1, umap_2, color = Neuron.type, shape = Neuron_type1), size = 1.5, alpha = 0.8) +
    scale_color_manual(values = CTCOLORS, name = "") +
    scale_shape_manual(values = c(21,23), name = "") + 
    theme_void() + 
    guides(color = guide_legend(override.aes = list(size=4)), shape = guide_legend(override.aes = list(size=4)))
  )
legend <- cowplot::get_legend(p)
png(here(plot_dir, "legend_subtypes_umap.png"), width = 300, height = 1800, res = 350)
grid::grid.newpage()
grid::grid.draw(legend)
dev.off()

# For heatmap
(p <- df_umap %>% 
  ggplot(aes(umap_1, umap_2)) +
  geom_point(aes(color = Neuron.type), size = 1.5, shape = 16, alpha = 0.5) +
  scale_color_manual(values = CTCOLORS, name = "") +
  theme_void() + 
  guides(colour = guide_legend(override.aes = list(size=4))))
legend <- cowplot::get_legend(p)
png(here(plot_dir, "legend_subtypes_heatmap.png"), width = 300, height = 1200, res = 320)
grid::grid.newpage()
grid::grid.draw(legend)
dev.off()


# ---- UMAP from VMRs for all methods ----
umapRegionalMeanMethyl <- function(method) {
  # `n_top`: number of top regions
  # `method`: either "vseq", "vseq_cr", "scbs", "smwd" or "scmet"  ('vseq' not available since no ranking available)
  df_umap <- fread(here(read_dir, paste0("metadata_umap_regional_methyl_", method, "_seed2010.txt.gz")))
  df_umap %>% 
    ggplot(aes(umap_1, umap_2)) +
    # geom_point(aes(color = Neuron.type), size = 1.5, shape = 16, alpha = 0.5) +
    # scale_color_manual(values = CTCOLORS) +
    geom_point(aes(fill = Neuron.type, shape = Neuron_type1), size = 1.5, alpha = 0.8, stroke = 0.3, color = 'grey16') +
    scale_fill_manual(values = CTCOLORS) +
    scale_shape_manual(values = c(21,23)) + 
    theme_classic() +
    theme(axis.text  = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.line  = element_blank())
  ggsave(here(plot_dir, paste0("umap_regional_methyl_", method, ".png")), width = 6, height = 5)
}

umapRegionalMeanMethyl(method = 'vseq')
umapRegionalMeanMethyl(method = 'vseq_cr')
umapRegionalMeanMethyl(method = 'scbs')
umapRegionalMeanMethyl(method = 'smwd')
umapRegionalMeanMethyl(method = 'scmet')

# ---- Heatmap of vmrseq top VMRs ----

res_ct <- list(
  'vseq' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_vseq.rds'))),
  'vseq_cr' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_vseq_cr.rds'))),
  'scbs' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_scbs.rds'))),
  'smwd' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_smwd.rds'))),
  'scmet' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_scmet.rds')))
)

# all(granges(res_region$vseq)$loglik_diff == granges(res_ct$vseq)$loglik_diff) # sanity check

heatmapTopRegions <- function(n_top, method, height, width, dissim_metric = "manhattan", hclust_metric = "ward") {
  # `n_top`: number of top regions
  # `method`: either "vseq", "scbs", "smwd" or "scmet"  ('vseq' not available since no ranking available)
  
  metric <- switch(method,
                   "vseq" = granges(res_region[[method]])$loglik_diff,
                   "scbs" = granges(res_region[[method]])$mcols.var,
                   "smwd" = granges(res_region[[method]])$var_lb,
                   "scmet" = granges(res_region[[method]])$tail_prob)

  ## Pick out top regions
  top_ind <- order(metric, decreasing = TRUE)[1:n_top]
  top.se <- res_region[[method]][top_ind]
  top_MF <- as.matrix(assays(top.se)$M / assays(top.se)$Cov)

  # Sort matrix and annotation in `Neuron.type` order
  idx <- order(md$Neuron.type)
  md_sorted <- md[idx, ]; rownames(md_sorted) <- md_sorted$sample
  top_MF_sorted <- top_MF[,idx]; colnames(top_MF_sorted) <- md_sorted$sample
  
  # Get hierarchical clustering on rows (i.e., regions)
  d_mat_rows <- as.matrix(cluster::daisy(top_MF_sorted, metric = dissim_metric, stand = FALSE))
  d_mat_rows[is.na(d_mat_rows)] <- mean(d_mat_rows, na.rm = TRUE)
  cluster_rows <- cluster::agnes(d_mat_rows, diss = TRUE, method = hclust_metric) %>% as.hclust()
  # ord <- cluster_rows$order
  # plot_MF <- top_MF_sorted[ord, ]
  # diff_exc_inh <- rowMeans(plot_MF[, md_sorted$Neuron_type1=='Excitatory'], na.rm = T) - rowMeans(plot_MF[, md_sorted$Neuron_type1=='Inhibitory'], na.rm = T)
  # clust_1 <- which(diff_exc_inh <= -0.1)
  # clust_2 <- which(diff_exc_inh >= 0.1)
  # clust_3 <- which(diff_exc_inh > -0.1 & diff_exc_inh < 0.1)
  # plot_MF <- plot_MF[c(clust_1, clust_2, clust_3),]
  
  graphics.off()
  png(here(plot_dir, paste0("heatmap_top", n_top, "regions_", method, ".png")), height = height, width = width, res = 170)
  # heatmap_palette <- colorRampPalette(RColorBrewer::brewer.pal(9, name = "GnBu")[1:7])(21)
  heatmap_palette <- colorRampPalette(RColorBrewer::brewer.pal(9, name = "GnBu")[c(4,5,6,8)])(21)
  # heatmap_palette <- colorRampPalette(viridis::inferno(9)[c(4:8)])(21)
  pheatmap::pheatmap(
    # plot_MF,
    top_MF_sorted,
    border_color = NA,
    color = heatmap_palette,
    cluster_cols = F,
    cluster_rows = cluster_rows,
    treeheight_row = 20,
    gaps_col = c(which(!duplicated(md_sorted$Neuron_type1)))-1,
    show_colnames = F,
    show_rownames = F,
    na_col = "white",
    # na_col = "grey90",
    annotation_col = md_sorted %>% select(Neuron.type, Neuron_type1),
    annotation_colors = list(Neuron.type = CTCOLORS, Neuron_type1 = BTCOLORS),
    annotation_legend = F,
    annotation_names_col = F
  )
  dev.off()
} 

heatmapTopRegions(n_top = 500, method = 'vseq', height = 800, width = 600)
heatmapTopRegions(n_top = 500, method = 'scbs', height = 800, width = 600)
heatmapTopRegions(n_top = 500, method = 'smwd', height = 800, width = 600)
heatmapTopRegions(n_top = 500, method = 'scmet', height = 800, width = 600)


# ---- Nearest neighbor scores ----

nnScorePlot <- function(k, theta, ylim, ybreaks) {
  
  score_broad.df <- fread(here(read_dir, paste0("nearest_neighbor_score_broadCellType_k", k, "_theta", theta, ".csv")))
  score_sub.df   <- fread(here(read_dir, paste0("nearest_neighbor_score_subCellType_k", k, "_theta", theta, ".csv")))
  score.df <- rbind(data.frame(score_broad.df, Label = 'Broad\nClasses'),
                    data.frame(score_sub.df,   Label = 'Subtypes'))
  p <- score.df %>%
    filter(Method != 'vmrseq CRs') %>%
    ggplot(aes(x = NTopRegions, y = NNScore, color = Method, shape = Label, linetype = Label)) +
    geom_point(size = 3) +
    geom_path() + 
    geom_point(data = score.df %>% filter(Method == 'vmrseq CRs'), 
               aes(x = NTopRegions, y = NNScore, color = Method, shape = Label), 
               size = 3) +
    scale_color_manual(values = COLORVALUES) +
    scale_shape_manual(values = c(1, 16)) +
    scale_linetype_manual(values = c(2, 1)) +
    scale_x_log10(breaks = c(300, 1000, 3000, 10000, 30000, 100000), labels = scales::comma) +
    scale_y_continuous(breaks = ybreaks, limits = ylim) +
    xlab("N Top Regions") + ylab("Nearest Neighbor Count Score") +
    guides(shape = guide_legend(title = 'Labeled by'),
           linetype = guide_legend(title = 'Labeled by')) +
    theme_classic()
  p
  ggsave(here(plot_dir, paste0("point_nnScore_vs_nTopRgions_k",k,"_theta",theta,".png")), height = 3.5, width = 4.5)
  return(p)
}

nnScorePlot(k = 100, theta = 0.7, ylim = c(0, 1), ybreaks = seq(0, 1, 0.5))

# ---- Percentage of VMLs out of all CpGs from homogeneous study ----

n_sites <- sapply(res_region, function(se) findOverlaps(GenomicRanges::reduce(granges(se)), sites.gr) %>% length())
pct <- n_sites / length(sites.gr)
ggplot() + 
  geom_bar(aes(x = methods, y = pct, color = methods, fill = methods), stat = 'identity') + 
  scale_y_continuous(labels = percent, limits = c(0, 0.118), name = '# VMLs / # CpGs in total (homogeneous population)') +
  scale_color_manual(values = COLORVALUES) + 
  scale_fill_manual(values = COLORVALUES) + 
  xlab('Methods') + 
  theme_classic() +
  theme(legend.position = 'none')
ggsave(here(plot_dir, 'barplot_pctSites_vs_methods_homStudy.png'), height = 4, width = 4)  

# # ---- Individual VMR ----
# 
# se <- loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs"))
# 
# plotCellTypeMeth <- function(top_i) {
#   
#   # index of the top ith region in se
#   idx <- order(granges(se)$loglik_diff, decreasing = T)[top_i]
#   gr <- granges(se)[idx]
#   
#   # summarize phenotypic info and methylation value per cell type
#   df <- data.frame(
#     md,
#     MF = assays(se)$M[idx,] / assays(se)$Cov[idx,]
#   ) %>% 
#     filter(!is.na(MF))
#   
#   # sub_title <- paste0("No.", top_i," Top Region")
#   main_title <- paste0(seqnames(gr), ": ", 
#                       format(start(gr), big.mark=","), "-", 
#                       format(end(gr), big.mark=","),
#                       " (width = ", format(end(gr)-start(gr), big.mark=","), " bp)")
#   df %>% ggplot(aes(Neuron.type, MF)) +
#     geom_jitter(height = 0.02, size = 0.05) +
#     # geom_jitter(aes(color = Neuron.type), height = 0.02, size = 0.05, alpha = 0.7) +
#     geom_violin(aes(fill = Neuron.type), alpha = 0.5, width = 1.3) +
#     scale_fill_manual(values = CTCOLORS) + 
#     labs(x = "Annotated Subtype",
#          y = "Regional Average Methylation",
#          fill = "Broad Neuron Type") +
#     facet_grid(~ Neuron_type1, scales = "free_x", space = "free_x") + 
#     ggtitle(main_title)  + 
#     guides(fill = "none") +
#     theme_classic() + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#           plot.title = element_text(hjust = 0.5),
#           plot.subtitle = element_text(hjust = 0.5))
#   ggsave(here(plot_dir, paste0("region_topNo", top_i, ".png")), width = 8, height = 3)
# } 
# 
# plotCellTypeMeth(top_i = 5)
# plotCellTypeMeth(top_i = 23)

# ---- cell type markers ----


plotCellTypeMeth <- function(method, idx, cell_type) {
  
  se <- res_region[[method]][idx,]
  gr <- granges(se)
  
  # summarize phenotypic info and methylation value per cell type
  df <- data.frame(
    md,
    MF = assays(se)$M / assays(se)$Cov
  ) %>% 
    filter(!is.na(MF))
  mean_meth_tg <- mean(df %>% filter(Neuron_type == cell_type) %>% pull(MF))
  mean_meth_bg <- mean(df %>% filter(Neuron_type != cell_type) %>% pull(MF))
  message('Average methylation in target cell type: ', round(mean_meth_tg, 2))
  message('Average methylation in background cell type: ', round(mean_meth_bg, 2))
  message('Difference in average methylation in traget and background cell type: ', round(mean_meth_bg - mean_meth_tg, 2))
  
  main_title <- methodName(method)
  sub_title <- paste0(seqnames(gr), ": ", 
                       format(start(gr), big.mark=","), "-", 
                       format(end(gr), big.mark=","),
                       " (width = ", format(end(gr)-start(gr)+1, big.mark=","), " bp)")
  p <- df %>% ggplot(aes(Neuron.type, MF)) +
    geom_jitter(height = 0.02, size = 0.05) +
    # geom_jitter(aes(color = Neuron.type), height = 0.02, size = 0.05, alpha = 0.7) +
    geom_violin(aes(fill = Neuron.type), alpha = 0.5, width = 1.3) +
    scale_fill_manual(values = CTCOLORS) + 
    scale_y_continuous(breaks = c(0,1)) +
    labs(x = "Annotated Subtype",
         y = "Regional Average\nMethylation",
         fill = "Broad Neuron Type") +
    facet_grid(~ Neuron_type1, scales = "free_x", space = "free_x") + 
    ggtitle(main_title, subtitle = sub_title)  + 
    guides(fill = "none") +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  return(p)
} 



plotTopCellTypeMarker <- function(method, cell_type, dmr_direction) {
  
  metric <- switch(method,
                   "vseq" = granges(res_region[[method]])$loglik_diff,
                   "scbs" = granges(res_region[[method]])$mcols.var,
                   "smwd" = granges(res_region[[method]])$var_lb,
                   "scmet" = granges(res_region[[method]])$tail_prob)
  
  # Makrer region indices
  mr_indices <- fread(here(read_dir, paste0('cellTypeMarkerRegionIndices_', method, '.csv'))) %>%
    filter(cellType == cell_type, dmrDirection == dmr_direction) %>%
    pull(dmrIndex)
  
  idx <- mr_indices[which.max(metric[mr_indices])]
  
  plotCellTypeMeth(method, idx, cell_type) 
  ggsave(here(plot_dir, paste0("violin_topRankedCellTypeMarker_", cell_type, "_", dmr_direction, "_", method,".png")), width = 8, height = 3)
}

plotTopCellTypeMarker(method = 'vseq', cell_type = 'mL4', dmr_direction = 'hypo')
# Average methylation in target cell type: 0.13
# Average methylation in background cell type: 0.75
# Difference in average methylation in traget and background cell type: 0.62
plotTopCellTypeMarker(method = 'scbs', cell_type = 'mL4', dmr_direction = 'hypo')
# Average methylation in target cell type: 0.1
# Average methylation in background cell type: 0.6
# Difference in average methylation in traget and background cell type: 0.5
plotTopCellTypeMarker(method = 'smwd', cell_type = 'mL4', dmr_direction = 'hypo')
# Average methylation in target cell type: 0.37
# Average methylation in background cell type: 0.73
# Difference in average methylation in traget and background cell type: 0.37
plotTopCellTypeMarker(method = 'scmet', cell_type = 'mL4', dmr_direction = 'hypo')
# Average methylation in target cell type: 0.55
# Average methylation in background cell type: 0.88
# Difference in average methylation in traget and background cell type: 0.33



# ---- variance across cell types of each method ----
ct.se <- list(
  'vseq' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_vseq.rds'))),
  'vseq_cr' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_vseq_cr.rds'))),
  'scbs' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_scbs.rds'))),
  'smwd' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_smwd.rds'))),
  'scmet' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_scmet.rds')))
)

var.list <- lapply(ct.se, function(SE) rowVars(as.matrix(assays(SE)$reginal_methyl), na.rm = TRUE))
var.df <- do.call(
  rbind, 
  lapply(1:length(var.list), function(i) data.frame(var = var.list[[i]], method = names(var.list)[i]))
) %>% 
  mutate(method = factor(method)) %>%
  mutate(method = fct_recode(method, 
                             "vmrseq" = "vseq", "vmrseq CRs" = "vseq_cr", 
                             "scbs" = "scbs", "Smallwood" = "smwd",
                             "scMET" = "scmet"))

methods <- c("vmrseq", "vmrseq CRs", "scbs", "Smallwood", "scMET")
var.df$method <- factor(var.df$method, levels = methods)

var.df %>%
  ggplot(aes(method, var, color = method)) + 
  geom_boxplot() +
  scale_color_manual(values = COLORVALUES) +
  ylim(0, 0.25) + 
  xlab("Method") + ylab("Variance of regional methyl across cell types") + 
  theme_classic()
ggsave(here(plot_dir, 'boxplot_cellTypeVariance_regionalMethylation.png'), width = 6, height = 4)

# ---- number of VMRs/VMLs detected by each methods ----


# Number of detected regions
n_regions <- sapply(res_region, function(se) granges(se) %>% GenomicRanges::reduce() %>% length())
names <- factor(sapply(names(n_regions), methodName), levels = sapply(names(n_regions), methodName))
ggplot() +
  geom_bar(aes(x = names, 
               y = n_regions, 
               color = names, 
               fill = names), stat = "identity") +
  scale_color_manual(values = COLORVALUES) +
  scale_fill_manual(values = COLORVALUES) +
  scale_y_continuous(labels = scales::scientific, name = "Number of regions detected") +
  xlab("Methods") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_line(color = "light grey"))
ggsave(here(plot_dir, "barplot_nNonOlapRgions_vs_methods.png"), height = 4, width = 5.5)

# Percentage of sites in detected regions
n_sites <- sapply(res_region, function(se) findOverlaps(GenomicRanges::reduce(granges(se)), sites.gr) %>% length())
pct <- n_sites / length(sites.gr)
ggplot() +
  geom_bar(aes(x = names, y = pct, color = names, fill = names), stat = "identity") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.1185), name = "Percentage of CpGs in detected regions") +
  scale_color_manual(values = COLORVALUES) +
  scale_fill_manual(values = COLORVALUES) +
  xlab("Methods") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_line(color = "light grey"))
ggsave(here(plot_dir, "barplot_pctSites_vs_methods.png"), height = 4, width = 4)


# Distribution of CpG number per detected region
n_sites_avail <- lapply(res_region, function(se) countOverlaps(granges(se), sites.gr))
ggplot() +
  geom_density(aes(n_sites_avail[[1]], color = names(COLORVALUES)[1])) +
  geom_density(aes(n_sites_avail[[2]], color = names(COLORVALUES)[2])) +
  geom_density(aes(n_sites_avail[[3]], color = names(COLORVALUES)[3])) +
  geom_density(aes(n_sites_avail[[4]], color = names(COLORVALUES)[4])) +
  geom_density(aes(n_sites_avail[[5]], color = names(COLORVALUES)[5])) +
  scale_x_continuous(limits = c(0, 400), name = "Number of CpG per detected region") +
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "light grey")) +
  scale_color_manual(name = "Methods", breaks = names(COLORVALUES), values = COLORVALUES)
ggsave(here(plot_dir, "densityplot_nCpGsPerRegion.png"), height = 3, width = 5)


# Distribution of width of detected region
region_width <- lapply(res_region, function(se) width(GenomicRanges::reduce(granges(se), min.gapwidth = 0)))
df <- map_df(region_width, 
             function(x) data.frame(Width = x), 
             .id = "Method") %>%
  mutate(Method = factor(sapply(Method, methodName), levels = methods))
df %>%
  ggplot(aes(Method, Width, color = Method)) +
  geom_boxplot() +
  scale_y_continuous(name = "Width of detected regions (bp)", limits = c(0, 40000)) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = "black"), 
        panel.grid = element_line(color = "light grey")) +
  scale_color_manual(name = "Methods", breaks = names(COLORVALUES), values = COLORVALUES)
ggsave(here(plot_dir, "densityplot_regionWidth.png"), height = 4, width = 4.5)

# Overlap between detected regions from methods
res.gr.list <- lapply(res_region, granges)
names(res.gr.list) <- sapply(names(res.gr.list), methodName)
comb_mat <- ComplexHeatmap::make_comb_mat(res.gr.list)
png(here(plot_dir, "UpSetPlot_allMethods_overlapWindowSize.png"), width = 1000, height = 600, res = 200)
ComplexHeatmap::UpSet(comb_mat, 
                      comb_order = order(ComplexHeatmap::comb_size(comb_mat), decreasing = T),
                      set_order = order(ComplexHeatmap::set_size(comb_mat), decreasing = F))
dev.off()

