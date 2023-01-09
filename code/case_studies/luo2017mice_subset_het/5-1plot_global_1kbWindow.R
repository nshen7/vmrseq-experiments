.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(scales)
library(pheatmap)

read_dir <- "data/interim/case_studies/luo2017mice_subset_het/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_subset_het/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_subset_het/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

md <- fread("data/interim/case_studies/luo2017mice_subset_het/metadata_luo2017mice_subset_het.csv")
rownames(md) <- md$sample

# ==== read in results from various methods ====
sites.gr <- readRDS(paste0(read_dir, "cpg_sites.rds"))
res_region <- list(
  vseq = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs_1kbWindow")),
  vseq_cr = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs_1kbWindow")),
  scbs = loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs")),
  smwd = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs")),
  scmet = loadHDF5SummarizedExperiment(paste0(read_dir, "scmet_regionSummary_vmrs"))
)
res_site <- list(
  vseq = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_siteSummary_sparseRep_vmrs_1kbWindow")),
  vseq_cr = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_siteSummary_sparseRep_crs_1kbWindow")),
  scbs = loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_siteSummary_sparseRep_vmrs")),
  smwd = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_siteSummary_sparseRep_vmrs")),
  scmet = loadHDF5SummarizedExperiment(paste0(read_dir, "scmet_siteSummary_sparseRep_vmrs"))
)

methods <- c("vmrseq", "CR in vmrseq", "scbs", "smallwood", "scMET")
methods <- factor(methods, levels = methods)

################################
###### General exploration #####
################################

# Color settings
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = COLORS[1], "CR in vmrseq" = COLORS[2], 
                  "scbs" = COLORS[3], "smallwood" = COLORS[4], "scMET" = COLORS[5])

# Number & percentage of sites in detected regions
n_sites <- sapply(res_region, function(se) findOverlaps(GenomicRanges::reduce(granges(se)), sites.gr) %>% length())
ggplot() +
  geom_bar(aes(x = methods, y = n_sites, color = methods, fill = methods), stat = "identity") +
  scale_y_continuous(labels = scientific, limits = c(0, 2.5e6), name = "Total number of CpGs in detected regions") +
  scale_color_manual(values = COLORVALUES) +
  scale_fill_manual(values = COLORVALUES) +
  xlab("Methods") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_line(color = "light grey"))
ggsave(paste0(plot_dir, "barplot_nSites_vs_methods_1kbWindow.png"), height = 4, width = 4)

# Percentage of sites in detected regions
pct <- n_sites / length(sites.gr)
ggplot() +
  geom_bar(aes(x = methods, y = pct, color = methods, fill = methods), stat = "identity") +
  scale_y_continuous(labels = percent, limits = c(0, 0.118), name = "Percentage of CpGs in detected regions") +
  scale_color_manual(values = COLORVALUES) +
  scale_fill_manual(values = COLORVALUES) +
  xlab("Methods") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_line(color = "light grey"))
ggsave(paste0(plot_dir, "barplot_pctSites_vs_methods_1kbWindow.png"), height = 4, width = 4)


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
ggsave(paste0(plot_dir, "densityplot_nCpGsPerRegion_1kbWindow.png"), height = 4, width = 7)


# Number of covered regions relative to total number of detected regions per cell
n_regn_avail <- lapply(res_region, function(se) colSums(assays(se)$Cov > 0)/nrow(se))
ggplot() +
  geom_density(aes(n_regn_avail[[1]], color = names(COLORVALUES)[1])) +
  geom_density(aes(n_regn_avail[[2]], color = names(COLORVALUES)[2])) +
  geom_density(aes(n_regn_avail[[3]], color = names(COLORVALUES)[3])) +
  geom_density(aes(n_regn_avail[[4]], color = names(COLORVALUES)[4])) +
  geom_density(aes(n_regn_avail[[5]], color = names(COLORVALUES)[5])) +
  scale_x_continuous(labels = percent, limits = c(0, 1), name = "Percentage of covered regions per cell") +
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid = element_line(color = "light grey")) +
  scale_color_manual(name = "Methods", breaks = names(COLORVALUES), values = COLORVALUES)
ggsave(paste0(plot_dir, "densityplot_pctCoveredRegionPerCell_1kbWindow.png"), height = 4, width = 7)


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
  
  # Sort matrix and annotation in `Neuron_type3` order
  idx <- order(md$Neuron_type3)
  md_sorted <- md[idx, ]; rownames(md_sorted) <- md_sorted$sample
  top_MF_sorted <- top_MF[,idx]; colnames(top_MF_sorted) <- md_sorted$sample
  
  # Get hierarchical clustering on rows (i.e., regions)
  d_mat_rows <- as.matrix(cluster::daisy(top_MF_sorted, metric = dissim_metric, stand = FALSE))
  d_mat_rows[is.na(d_mat_rows)] <- mean(d_mat_rows, na.rm = TRUE)
  cluster_rows <- cluster::agnes(d_mat_rows, diss = TRUE, method = hclust_metric)
  
  graphics.off() 
  png(paste0(plot_dir, "heatmap_top", n_top, "regions_", method, "_1kbWindow.png"), height = 900, width = 900)
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
           annotation_col = md_sorted %>% select(Neuron_type3, Neuron_type1),
           main = paste0("Regional Mean Methylation of Top ", n_top, " Regions from ", method_name)
  )
  dev.off()
} 

heatmapTopRegions(n_top = 100, method = "vseq")
heatmapTopRegions(n_top = 300, method = "vseq")
heatmapTopRegions(n_top = 1000, method = "vseq")

###################################################################
###### Clustering analysis based on regional mean methylation #####
###################################################################

umapRegionalMeanMethyl <- function(method, dissim_metric = "manhattan", umap_metric = "euclidean") {
  # `n_top`: number of top regions
  # `method`: either "vseq", "vseq_cr", "scbs", "smwd" or "scmet"  ('vseq' not available since no ranking available)
  method_name <- switch(method,
                        "vseq" = "vmrseq",
                        "vseq_cr" = "vmrseq CRs",
                        "scbs" = "scbs",
                        "smwd" = "Smallwood",
                        "scmet" = "scMET")
  set.seed(2022)
  
  MF <- as.matrix(assays(res_region[[method]])$M/assays(res_region[[method]])$Cov) 

  # Get dissimilarity matrix for columns (i.e., cells)
  d_mat_cols <- cluster::daisy(t(MF), metric = dissim_metric, stand = FALSE) %>% as.matrix()
  rownames(d_mat_cols) <- colnames(d_mat_cols) <- md$sample
  fwrite(d_mat_cols, paste0(write_dir, "dissimilarity_matrix_regional_methyl_", method, "_1kbWindow.txt.gz"),
         col.names = TRUE, row.names = TRUE, quote = F)
  
  # Get UMAP embedding
  umap <- uwot::umap(d_mat_cols %>% as("sparseMatrix"), n_neighbors = 15, n_components = 2, metric = umap_metric)
  df_umap <- data.frame(md, umap_1 = umap[,1], umap_2 = umap[,2])
  fwrite(df_umap, paste0(write_dir, "metadata_umap_regional_methyl_", method, "_1kbWindow.txt.gz"),
         col.names = TRUE, row.names = FALSE, quote = F)
  
  df_umap %>% 
    ggplot(aes(umap_1, umap_2, color = Neuron_type3, shape = Neuron_type1)) +
    geom_point() +
    scale_color_brewer(palette="Set2") + 
    theme_classic() +
    ggtitle(paste0("Method: ", method_name))
  ggsave(paste0(plot_dir, "umap_regional_methyl_", method, "_1kbWindow.png"), width = 7, height = 5.5)
}

umapRegionalMeanMethyl("vseq")
umapRegionalMeanMethyl("vseq_cr")


#################################################################
###### Clustering analysis based on single-site methylation #####
#################################################################

umapSiteMethyl <- function(method, dissim_metric = "manhattan", umap_metric = "euclidean") {
  # `n_top`: number of top regions
  # `method`: either "vseq", "scbs", "smwd" or "scmet"  ('vseq' not available since no ranking available)
  method_name <- switch(method,
                        "vseq" = "vmrseq",
                        "vseq_cr" = "vmrseq CRs",
                        "scbs" = "scbs",
                        "smwd" = "Smallwood",
                        "scmet" = "scMET")
  set.seed(2022)
  
  M_mat <- assays(res_site[[method]])$M_mat %>% as("sparseMatrix") %>% recommenderlab::dropNA2matrix()
  
  # Get UMAP embedding for columns (i.e., cells)
  d_mat_cols <- cluster::daisy(t(M_mat), metric = dissim_metric, stand = FALSE)%>% as.matrix()
  rownames(d_mat_cols) <- colnames(d_mat_cols) <- md$sample
  fwrite(d_mat_cols, paste0(write_dir, "dissimilarity_matrix_site_methyl_", method, "_1kbWindow.txt.gz"),
         col.names = TRUE, row.names = TRUE, quote = F)
  
  # Get UMAP embedding
  umap <- uwot::umap(d_mat_cols %>% as("sparseMatrix"), n_neighbors = 15, n_components = 2, metric = umap_metric)
  df_umap <- data.frame(md, umap_1 = umap[,1], umap_2 = umap[,2])
  fwrite(df_umap, paste0(write_dir, "metadata_umap_site_methyl_", method, "_1kbWindow.txt.gz"),
         col.names = TRUE, row.names = FALSE, quote = F)
  
  df_umap %>% 
    ggplot(aes(umap_1, umap_2, color = Neuron_type3, shape = Neuron_type1)) +
    geom_point() +
    scale_color_brewer(palette="Set2") + 
    theme_classic() +
    ggtitle(paste0("Method: ", method_name))
  ggsave(paste0(plot_dir, "umap_site_methyl_", method, "_1kbWindow.png"), width = 7, height = 5.5)
}

umapSiteMethyl("vseq")
umapSiteMethyl("vseq_cr")

