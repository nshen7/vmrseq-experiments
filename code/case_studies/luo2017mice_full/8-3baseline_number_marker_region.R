source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(pheatmap)

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary"
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

# ---- Fixed arguments ----
min_diff_hyper <- 0.2; min_diff_hypo <- 0.2
seeds <- 1:10

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

# ---- Load data ----
md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")
res.se <- list(
  'vseq' = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs")),
  'vseq_cr' = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs")),
  'scbs' = loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs")),
  'smwd' = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs")),
  'scmet' = loadHDF5SummarizedExperiment(paste0(read_dir, "scmet_regionSummary_vmrs"))
)

# ---- Summarize # marker regions from shuffled cell type ----
countMarkerRegionsShuffleCellType <- function(target_type, method, min_diff_hyper, min_diff_hypo, min_avail_cell_prop = 0.2) {
  
  myCount <- function(seed) {
    ct_method.se <- readRDS(here(read_dir, 'cell_type_shuffled',
                                 paste0('SummarizedExperiment_shuffledCellType_regionalMean_', method, '_seed', seed, '.rds')))
    list_idx <- findCellTypeMarkerIndices(
      target_type = target_type, 
      target_type_n_cell = sum(md$Neuron.type == target_type), 
      ct.se = ct_method.se, 
      min_diff_hyper = min_diff_hyper, 
      min_diff_hypo = min_diff_hypo, 
      min_avail_cell_prop = min_avail_cell_prop
    )
    df <- data.frame(method     = method,
                     cellType  = target_type,
                     markerProp = c(length(list_idx$idx_hyper), length(list_idx$idx_hypo))/nrow(res.se[[method]]),
                     dmrDirection  = c('hyper', 'hypo'),
                     seed       = seed)
    return(df)
  }
  
  counts.df <- map_dfr(seeds, myCount)
  return(counts.df)
}

for (method in names(res.se)) {
  counts_method.df <- map_dfr(.x = unique(md$Neuron.type), .f = countMarkerRegionsShuffleCellType,
                            method = method,
                            min_diff_hyper = min_diff_hyper,
                            min_diff_hypo = min_diff_hypo)
  fwrite(counts_method.df, here(write_dir, paste0('markerRegionCounts_shuffledCellTypeLabels_', method, '.csv')))
  message('Finished ', method, '...')
}


# ---- Plot fake n_marker from shuffling cell type labels and real n_marker ----

counts_real_all.df <- fread(here(write_dir, paste0('markerRegionCounts_realCellTypeLabels_allMethods.csv')))

for (this_method in names(res.se)) {
  counts_real.df <- counts_real_all.df  %>%
    arrange(cellProp) %>%  # order cell types by cellProp
    mutate(cellType = factor(cellType, levels = unique(cellType))) %>% 
    filter(method == methodName(this_method))
  counts_fake.df <- fread(here(write_dir, paste0('markerRegionCounts_shuffledCellTypeLabels_', this_method, '.csv'))) %>%
    mutate(cellType = factor(cellType, levels = levels(counts_real.df$cellType)))
  
  for (direction in c('hyper', 'hypo')) {
    scale <- max(counts_real.df$cellProp) / max(counts_real.df$markerProp)
    counts_real_sub.df <- counts_real.df %>%
      filter(dmrDirection == direction & method == methodName(this_method))
    ggplot() +
      geom_boxplot(data = counts_fake.df, aes(cellType, markerProp), color = 'grey') +
      geom_line(data = counts_real_sub.df, aes(cellType, markerProp, color = method, group = method)) +
      scale_color_manual(values = COLORVALUES) +
      geom_point(data = counts_real_sub.df, aes(cellType, cellProp/scale, shape = 'Cell proportion'), color = 'darkgreen') +
      scale_shape_manual(values =  c('Cell proportion' = 4)) +
      scale_y_continuous(sec.axis = sec_axis(~ .*scale, name = 'Cell proportion')) +
      theme_classic() +
      ylab('# markers / # VMRs') + xlab('Cell type')
    ggsave(here(plot_dir, paste0('markerProp_vs_cellType_comparedToShuffledLabel_orderedByCellProp_', direction, '_', this_method, '.png')), width = 10, height = 6)
  }
  message('Finished ', this_method, '...')
}

