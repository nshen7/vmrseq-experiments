source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(pheatmap)
source("code/case_studies/util_functions/8findCellTypeMarkers.R")

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

# ---- fixed arguments ----
min_diff_hyper <- 0.2; min_diff_hypo <- 0.2

# ---- load data ----
res.se <- list(
  'vseq' = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs")),
  'vseq_cr' = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs")),
  'scbs' = loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs")),
  'smwd' = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs")),
  'scmet' = loadHDF5SummarizedExperiment(paste0(read_dir, "scmet_regionSummary_vmrs"))
)


# Read in SE object of regional summary of VMRs 
md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")

# Color settings
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = COLORS[1], "vmrseq CRs" = COLORS[2],
                 "scbs" = COLORS[3], "Smallwood" = COLORS[4], "scMET" = COLORS[5])
methodName <- function(method) switch (method,
                                       'vseq' = 'vmrseq',
                                       'vseq_cr' = 'vmrseq CRs',
                                       'scbs' = 'scbs',
                                       'smwd' = 'Smallwood',
                                       'scmet' = 'scMET')


# ---- Compare distribution of variance across cell types ----

## Regional average methylation for each cell type
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
ggsave(here(plot_dir, 'boxplot_cellTypeVariance_regionalMethylation.png'), width = 6, height = 5)

var.df %>%
  mutate(var = ifelse(var > 0.15, 0.15, var)) %>%
  ggplot(aes(x = var, y = ..density.., fill = method)) + 
  geom_histogram(color = "white", binwidth = 0.002, boundary = 0) +
  facet_wrap(~ method, ncol = 1) + 
  scale_fill_manual(values = COLORVALUES) +
  xlim(0, 0.152) + ylim(0, 125) + 
  theme_classic()
ggsave(here(plot_dir, 'histogram_cellTypeVariance_regionalMethylation.png'), width = 4, height = 6)

var.df %>%
  filter(var < 0.02) %>%
  ggplot(aes(x = var, y = ..density.., fill = method)) + 
  geom_histogram(color = "white") +
  facet_wrap(~ method) + 
  scale_fill_manual(values = COLORVALUES) +
  xlim(0, 0.02) + 
  theme_classic()
ggsave(here(plot_dir, 'histogram_cellTypeVariance_regionalMethylation_leq0.02.png'), width = 7, height = 5)


# ---- Identify cell-type-specific marker regions from VMRs ----
ct.se <- list(
  'vseq' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_vseq.rds'))),
  'vseq_cr' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_vseq_cr.rds'))),
  'scbs' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_scbs.rds'))),
  'smwd' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_smwd.rds'))),
  'scmet' = readRDS(here(read_dir, paste0('SummarizedExperiment_cellType_regionalMean_scmet.rds')))
)

findMarkerRegions <- function(target_type, method, min_diff_hyper, min_diff_hypo, min_avail_cell_prop = 0.2) {
  
  ct_method.se <- ct.se[[method]]
  
  list_idx <- findCellTypeMarkerIndices(
    target_type = target_type, 
    target_type_n_cell = sum(md$Neuron.type == target_type), 
    ct.se = ct_method.se, 
    min_diff_hyper = min_diff_hyper, 
    min_diff_hypo = min_diff_hypo, 
    min_avail_cell_prop = min_avail_cell_prop
  )
  idx_hyper <- list_idx$idx_hyper; idx_hypo <- list_idx$idx_hypo

  if (length(idx_hyper) > 0) {
    idx_hyper.df <- data.frame(dmrIndex = idx_hyper,
                               dmrDirection  = rep('hyper', length(idx_hyper)),
                               cellType  = target_type)
  } else {
    idx_hyper.df <- NULL
  }
  
  if (length(idx_hypo) > 0) {
    idx_hypo.df <- data.frame(dmrIndex = idx_hypo,
                              dmrDirection  = rep('hypo', length(idx_hypo)),
                              cellType  = target_type)
  } else {
    idx_hypo.df <- NULL
  }
  
  if (length(idx_hyper) > 0 | length(idx_hypo) > 0) {
    return(rbind(idx_hyper.df, idx_hypo.df))
  } else {
    return(NULL)
  }
}

for (method in names(ct.se)) {
  idx_all.df <- map_dfr(.x = colnames(ct.se[[method]]), .f = findMarkerRegions, 
                        method = method, 
                        min_diff_hyper = min_diff_hyper, min_diff_hypo = min_diff_hypo)
  fwrite(idx_all.df, here(write_dir, paste0('cellTypeMarkerRegionIndices_', method, '.csv')))
}

# ---- Plot heatmap of hypermethylated cell-type marker regions ----

plotMarkerRegions <- function(method, direction) {
  
  ## Add metadata into SE object
  method.se <- res.se[[method]]
  
  # all(colnames(method.se) == md$sample) # = TRUE
  colData(method.se) <- cbind(colData(method.se), md)
  
  idx_all.df <- fread(here(write_dir, paste0('cellTypeMarkerRegionIndices_', method, '.csv')))
  
  # Sort SE in `Neuron.type` order
  cell_order <- colData(method.se) %>%
    as.data.frame() %>%
    mutate(order = 1:ncol(method.se)) %>%
    arrange(Neuron_type1, Neuron.type) %>%
    pull(order)
  
  # Sort rows of cell-type-specific VMRs in `Neuron.type` order
  marker_idx.df <- idx_all.df %>%
    filter(dmrDirection == direction, ) %>%
    arrange(match(cellType, unique(colData(method.se)$Neuron.type[cell_order])))
  
  # Extract marker regions with cells in Neuron.type order
  marker.se <- method.se[marker_idx.df$dmrIndex, cell_order]
  
  # Sort cols of cells in `Neuron.type` order
  marker_MF <- as.matrix(assays(marker.se)$M / assays(marker.se)$Cov)
  colnames(marker_MF) <- colData(marker.se)$sample
  rownames(marker_MF) <- marker_idx.df$dmrIndex
  
  # Set palette for heatmap blocks
  heatmap_palette <- colorRampPalette(RColorBrewer::brewer.pal(4, name = 'PuOr'))(21)
  
  # Set palette for heatmap annotations
  ct <- unique(colData(marker.se)$Neuron.type)
  ann_palette <- colorRampPalette(RColorBrewer::brewer.pal(8, name = 'Set3'))(length(ct))
  names(ann_palette) <- ct
  ann_colors = list(
    Neuron.type = ann_palette
  )
  
  graphics.off()
  png(paste0(plot_dir, paste0('heatmap_cellType_', direction, 'Regions_', method, '.png')), height = 900, width = 900)
  pheatmap(marker_MF,
           border_color = NA,
           color = heatmap_palette,
           cluster_cols = F,
           cluster_rows = F,
           show_colnames = F,
           show_rownames = F,
           # na_col = 'grey90',
           na_col = 'white',
           annotation_col = colData(marker.se) %>% 
             as.data.frame() %>% 
             column_to_rownames('sample') %>% 
             select(Neuron.type, Neuron_type1),
           annotation_row = marker_idx.df %>%
             as.data.frame() %>% 
             column_to_rownames('dmrIndex') %>% 
             # select(cellType),
             mutate(Neuron.type = cellType) %>%
             select(Neuron.type),
           annotation_colors = ann_colors,
           main = paste0('Cell-type ', direction, '-methylated marker regions from ', method)
  )
  dev.off()
  message('Finished ', method, ' ', direction, '...')
}


for (method in names(res.se)) {
  for (direction in c('hyper', 'hypo')) {
    plotMarkerRegions(method, direction)
  }
}


# ---- Count marker regions per cell type per method ----

ct_prop_0.df <- md %>%
  group_by(Neuron.type) %>%
  summarise(cellProp = n()/nrow(md)) %>%
  dplyr::rename(cellType = Neuron.type)

ct_prop.df <- rbind(
  data.frame(dmrDirection = 'hyper', ct_prop_0.df),
  data.frame(dmrDirection = 'hypo', ct_prop_0.df)
)

countMarkerRegions <- function(this_method)  {
  marker_count.df <- fread(here(write_dir, paste0('cellTypeMarkerRegionIndices_', this_method, '.csv'))) %>%
    group_by(dmrDirection, cellType) %>% 
    summarise(markerProp = n()/nrow(res.se[[this_method]])) %>% 
    full_join(ct_prop.df, by = c('dmrDirection', 'cellType')) %>%
    mutate(method = methodName(this_method), markerProp = ifelse(is.na(markerProp), 0, markerProp))
  return(marker_count.df)
}

marker_count_all.df <- map_dfr(names(res.se), countMarkerRegions)
fwrite(marker_count_all.df, here(write_dir, paste0('markerRegionCounts_realCellTypeLabels_allMethods.csv')))

# Plot number of markers vs. cell proportion
for (direction in c('hyper', 'hypo')) {
  scale <- max(marker_count_all.df$cellProp) / max(marker_count_all.df$markerProp)
  marker_count_all.df %>%
    filter(dmrDirection == direction) %>%
    arrange(cellProp) %>%
    mutate(cellType = factor(cellType, levels = unique(cellType))) %>% # order cell types by cellProp
    ggplot(aes(cellType, markerProp, color = method, group = method)) +
    # geom_point(aes(shape = '# markers / # VMRs')) +
    geom_line() +
    geom_point(aes(cellType, cellProp/scale, shape = 'Cell proportion'), color = 'darkgreen') +
    # geom_line(aes(cellType, cellProp/scale, group = 'Cell proportion'), color = 'darkgreen') +
    scale_color_manual(values = COLORVALUES) +
    scale_y_continuous(sec.axis = sec_axis(~ .*scale, name = 'Cell proportion')) +
    scale_shape_manual(values =  c('Cell proportion' = 4)) +
    theme_classic() +
    ylab('# markers / # VMRs') + xlab('Cell type')
  ggsave(here(plot_dir, paste0('markerProp_vs_cellType_orderedByCellProp_', direction, '.png')), width = 10, height = 6)
}

