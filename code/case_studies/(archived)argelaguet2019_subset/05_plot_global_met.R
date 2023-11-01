source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(scales)
library(pheatmap)

read_dir_met <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', '04_summarize_output_met')
write_dir <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', '05_plot_global_met')
if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- here('plots', 'case_studies', 'argelaguet2019_subset', '05_plot_global_met')
if (!file.exists(plot_dir)) dir.create(plot_dir)

# ---- Load data ----
md <- fread(here('data', 'metadata', 'argelaguet2019', 'argelaguet2019_subset_met&rna_sample_metadata_processed.csv'))

res_region <- list(
  vseq = loadHDF5SummarizedExperiment(here(read_dir_met, "vmrseq_regionSummary_vmrs"))#,
  # vseq_cr = loadHDF5SummarizedExperiment(here(read_dir_met, "vmrseq_regionSummary_crs")),
  # scbs = loadHDF5SummarizedExperiment(paste0(read_dir_met, "scbs_regionSummary_vmrs")),
  # smwd = loadHDF5SummarizedExperiment(paste0(read_dir_met, "smallwood_regionSummary_vmrs")),
  # scmet = loadHDF5SummarizedExperiment(paste0(read_dir_met, "scmet_regionSummary_vmrs"))
)
dim(res_region$vseq) # 31154   726

###############################################################
###### Heatmap of methylation in highly variable regions ######
###############################################################

heatmapTopRegions <- function(n_top, method, dissim_metric = "manhattan", hclust_metric = "ward") {
  # `n_top`: number of top regions
  # `method`: either "vseq", "scbs", "smwd" or "scmet"  ('vseq' not available since no ranking available)
  
  # Extact top VMRs 
  metric <- switch(method,
                   "vseq" = granges(res_region[[method]])$loglik_diff,
                   "scbs" = granges(res_region[[method]])$mcols.var,
                   "smwd" = granges(res_region[[method]])$var_lb,
                   "scmet" = granges(res_region[[method]])$tail_prob)
  method_name <- switch(method,
                        "vseq" = "vmrseq",
                        "scbs" = "scbs",
                        "smwd" = "Smallwood",
                        "scmet" = "scMET")
  
  top_ind <- order(metric, decreasing = TRUE)[1:n_top]
  top.se <- res_region[[method]][top_ind]
  top_MF <- as.matrix(assays(top.se)$M / assays(top.se)$Cov)
  
  # Sort matrix and annotation in `lineage10x_2` order
  idx <- order(md$stage, md$lineage10x_2)
  md_sorted <- md[idx, ]; rownames(md_sorted) <- md_sorted$sample
  top_MF_sorted <- top_MF[,idx]; colnames(top_MF_sorted) <- md_sorted$sample
  
  # Get hierarchical clustering on rows (i.e., regions)
  d_mat_rows <- as.matrix(cluster::daisy(top_MF_sorted, metric = dissim_metric, stand = FALSE))
  d_mat_rows[is.na(d_mat_rows)] <- mean(d_mat_rows, na.rm = TRUE)
  cluster_rows <- cluster::agnes(d_mat_rows, diss = TRUE, method = hclust_metric)
  
  graphics.off() 
  png(here(plot_dir, paste0("heatmap_top", n_top, "regions_", method, ".png")), height = 900, width = 900)
  heatmap_palette <- colorRampPalette(RColorBrewer::brewer.pal(8, name = "YlOrRd"))(21)
  pheatmap(top_MF_sorted[cluster_rows$order, ],
           # cluster_rows = cluster_rows,
           # cluster_cols = cluster_results$hclust_obj,
           # treeheight_row = 0,
           border_color = NA,
           color = heatmap_palette,
           cluster_cols = F,
           cluster_rows = F,
           show_colnames = F,
           show_rownames = F,
           na_col = "grey90",
           annotation_col = md_sorted %>% select(lineage10x_2, stage),
           main = paste0("Regional Mean Methylation of Top ", n_top, " Regions from ", method_name)
  )
  dev.off()
} 

heatmapTopRegions(n_top = 1000, method = "vseq")
heatmapTopRegions(n_top = 3000, method = "vseq")

##################################################################
###### Clustering analysis based on regional mean methylation #####
###################################################################

# ---- Compute dissimilarity matrix ----
computeDissimMatrix <- function(method, 
                                top_n_regions = NULL,
                                dissim_metric = "manhattan", 
                                umap_metric = "euclidean",
                                seed = 2010) {
  # `n_top`: number of top regions
  # `method`: either "vseq", "vseq_cr", "scbs", "smwd" or "scmet"  ('vseq' not available since no ranking available)
  method_name <- switch(method,
                        "vseq" = "vmrseq",
                        "vseq_cr" = "vmrseq CRs",
                        "scbs" = "scbs",
                        "smwd" = "Smallwood",
                        "scmet" = "scMET")
  metric <- switch(method,
                   "vseq" = granges(res_region[[method]])$loglik_diff,
                   "scbs" = granges(res_region[[method]])$mcols.var,
                   "smwd" = granges(res_region[[method]])$var_lb,
                   "scmet" = granges(res_region[[method]])$tail_prob)
  if (!is.null(top_n_regions) & method=="vseq_cr") stop("CRs from vmrseq does not have rank.")
  
  set.seed(seed)
  
  se <- res_region[[method]]
  
  if (!is.null(top_n_regions)) { # take top n regions for clustering
    top_ind <- order(metric, decreasing = TRUE)[1:top_n_regions]
    se <- se[top_ind]
  }
  
  MF <- as.matrix(assays(se)$M/assays(se)$Cov) 
  
  name_seg <- ifelse(is.null(top_n_regions), yes = "", no = paste0("_top", top_n_regions, "regions"))
  
  # Get dissimilarity matrix for columns (i.e., cells)
  d_mat_cols <- cluster::daisy(t(MF), metric = dissim_metric, stand = FALSE) %>% as.matrix()
  rownames(d_mat_cols) <- colnames(d_mat_cols) <- md$sample
  fwrite(d_mat_cols, here(write_dir, paste0("dissimilarity_matrix_regional_methyl_", method, name_seg, ".txt.gz")),
         col.names = TRUE, row.names = TRUE, quote = F)
  
  # Get UMAP embedding
  umap <- uwot::umap(d_mat_cols %>% as("sparseMatrix"), n_neighbors = 15, n_components = 2, metric = umap_metric)
  df_umap <- data.frame(md, umap_1 = umap[,1], umap_2 = umap[,2])
  fwrite(df_umap, here(write_dir, paste0("metadata_umap_regional_methyl_", method, name_seg, "_seed", seed, ".txt.gz")),
         col.names = TRUE, row.names = FALSE, quote = F)
}
## Run
computeDissimMatrix(method = 'vseq', top_n_regions = 3000)
computeDissimMatrix(method = 'vseq')


# ---- Plot UMAP projection ----
umapRegionalMeanMethyl <- function(method,
                                   top_n_regions = NULL, 
                                   seed = 2010) {
  method_name <- switch(method,
                        "vseq" = "vmrseq",
                        "vseq_cr" = "vmrseq CRs",
                        "scbs" = "scbs",
                        "smwd" = "Smallwood",
                        "scmet" = "scMET")
  name_seg <- ifelse(is.null(top_n_regions), yes = "", no = paste0("_top", top_n_regions, "regions"))
  df_umap <- fread(here(write_dir, paste0("metadata_umap_regional_methyl_", method, name_seg, "_seed", seed, ".txt.gz")))
  
  STAGE_SHAPE <- c(4,15,16,17)
  names(STAGE_SHAPE) <- unique(df_umap$stage)
  
  PHASE_COLOR <- c("#0072B2", "#ABDDA4", "#CC79A7")
  names(PHASE_COLOR) <- levels(factor(unique(df_umap$Phase)))
  
  df_umap %>% 
    ggplot(aes(S.Score, G2M.Score, color = Phase)) + 
    geom_point() +
    scale_colour_manual(values = PHASE_COLOR) +
    theme_classic()
  ggsave(here(plot_dir, paste0("point_cellCycleScore_G2M_vs_S.png")),  width = 6, height = 5)
  
  ## Plot UMAP with stage, cell type and cell cycle phase colored
  df_umap %>% 
    ggplot(aes(umap_1, umap_2, color = stage, shape = stage)) +
    geom_point(size = 0.5) +
    scale_color_brewer(palette = "Paired") +
    scale_shape_manual(values = STAGE_SHAPE) + 
    theme_classic() + 
    ggtitle(paste0("Method: ", method_name))
  ggsave(here(plot_dir, paste0("umap_regional_methyl_coloredByStage_", method, name_seg, "_seed", seed, ".png")), width = 4, height = 3.5)
  
  for (this_stage in unique(df_umap$stage)) {
    df_umap %>% 
      filter(stage == this_stage) %>%
      ggplot(aes(umap_1, umap_2, color = lineage10x_2)) +
      geom_point(size = 0.3) +
      xlim(min(df_umap$umap_1), max(df_umap$umap_1)) + 
      ylim(min(df_umap$umap_2), max(df_umap$umap_2)) + 
      theme_classic() +
      ggtitle(paste0('Stage ', this_stage))
    ggsave(here(plot_dir, paste0("umap_regional_methyl_lineage10x_2_", method, name_seg, "_seed", seed, "_stage", this_stage, ".png")), width = 4, height = 3.5)
    df_umap %>% 
      filter(stage == this_stage) %>%
      ggplot(aes(umap_1, umap_2, color = lineage10x)) +
      geom_point(size = 0.3) +
      xlim(min(df_umap$umap_1), max(df_umap$umap_1)) + 
      ylim(min(df_umap$umap_2), max(df_umap$umap_2)) + 
      theme_classic() +
      ggtitle(paste0('Stage ', this_stage))
    ggsave(here(plot_dir, paste0("umap_regional_methyl_lineage10x_", method, name_seg, "_seed", seed, "_stage", this_stage, ".png")), width = 4, height = 3.5)
    df_umap %>% 
      filter(stage == this_stage) %>%
      ggplot(aes(umap_1, umap_2, color = Phase)) +
      geom_point(size = 0.3) +
      xlim(min(df_umap$umap_1), max(df_umap$umap_1)) + 
      ylim(min(df_umap$umap_2), max(df_umap$umap_2)) + 
      theme_classic() +
      scale_colour_manual(values = PHASE_COLOR) + 
      ggtitle(paste0('Stage ', this_stage))
    ggsave(here(plot_dir, paste0("umap_regional_methyl_cellCyclePhase_", method, name_seg, "_seed", seed, "_stage", this_stage, ".png")), width = 4, height = 3.5)
  }
  
  df_umap %>% 
    ggplot(aes(umap_1, umap_2, color = Phase)) +
    geom_point(size = 1) +
    scale_colour_manual(values = PHASE_COLOR) + 
    theme_classic() +
    ggtitle(paste0("Method: ", method_name))
  ggsave(here(plot_dir, paste0("umap_regional_methyl_cellCyclePhase_", method, name_seg, "_seed", seed, ".png")), width = 4, height = 3.5)
  
  for (this_phase in unique(df_umap$Phase)) {
    df_umap %>%
      filter(Phase == this_phase) %>%
      ggplot(aes(umap_1, umap_2)) +
      geom_point(size = 0.4, color = PHASE_COLOR[this_phase]) +
      xlim(min(df_umap$umap_1), max(df_umap$umap_1)) +
      ylim(min(df_umap$umap_2), max(df_umap$umap_2)) +
      theme_classic() +
      ggtitle(paste0('Phase ', this_phase))
    ggsave(here(plot_dir, paste0("umap_regional_methyl_cellCyclePhase_", method, name_seg, "_seed", seed, "_phase", this_phase, ".png")), width = 4, height = 3.5)
  }
  
  df_umap %>% 
    ggplot(aes(umap_1, umap_2, color = S.Score)) +
    geom_point(size = 1.3) +
    theme_classic() +
    scale_color_gradient(low = PHASE_COLOR['G1'], high = PHASE_COLOR['S']) +
    ggtitle(paste0("Method: ", method_name))
  ggsave(here(plot_dir, paste0("umap_regional_methyl_cellCycleScore_S_", method, name_seg, "_seed", seed, ".png")), width = 4, height = 3.5)
  df_umap %>% 
    ggplot(aes(umap_1, umap_2, color = G2M.Score)) +
    geom_point(size = 1.3) +
    scale_color_gradient(low = PHASE_COLOR['G1'], high = PHASE_COLOR['G2M']) +
    theme_classic() +
    ggtitle(paste0("Method: ", method_name))
  ggsave(here(plot_dir, paste0("umap_regional_methyl_cellCycleScore_G2M_", method, name_seg, "_seed", seed, ".png")), width = 4, height = 3.5)
}
## Run
umapRegionalMeanMethyl(method = 'vseq')
umapRegionalMeanMethyl(method = 'vseq', top_n_regions = 3000)
