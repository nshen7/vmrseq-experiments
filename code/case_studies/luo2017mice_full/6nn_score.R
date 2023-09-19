source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(here)
source("code/case_studies/util_functions/6nnScore.R")

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)
md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")

# ---- utils ----

methodName <- function(method) switch (method,
                                       'vseq' = 'vmrseq',
                                       'vseq_cr' = 'vmrseq CRs',
                                       'scbs' = 'scbs',
                                       'smwd' = 'Smallwood',
                                       'scmet' = 'scMET')
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = COLORS[1], "vmrseq CRs" = COLORS[2],
                 "scbs" = COLORS[3], "Smallwood" = COLORS[4], "scMET" = COLORS[5])

# ---- main ----

computeScoreForAll <- function(k, theta, cell_type) {
  
  
  if (cell_type == 'broad') {
    # Use main neuron types as reference cell labels: 'Excitatory' or 'inhibitory'
    true_clust <- md$Neuron_type1
  } else if (cell_type == 'sub') {
    # Use neuron subtypes as reference cell labels: 'mL5-1' etc
    true_clust <- md$Neuron.type
  }
  
  ## Nearest neighbor score computed from regional mean methyl for each method
  nnScoreMethod <- function(method, top_n_regions = NULL) {
    
    name_seg <- ifelse(top_n_regions == '', yes = '', no = paste0("_top", top_n_regions, "regions"))
    path <- paste0(write_dir, "dissimilarity_matrix_regional_methyl_", method, name_seg, ".txt.gz")
    if (file.exists(path)) d_mat <- fread(path, drop = 1) else stop('Dissimilarity matrix not written yet!')
    
    stopifnot(all(md$sample == names(d_mat)))
    return(nnScore(d_mat, true_clust, k, theta))
  }
  
  ## Summarize nn score from various top n regions
  score.df <- expand.grid(
    Method = c('vseq', 'scbs', 'smwd', 'scmet'), 
    NTopRegions = c('300', '1000', '3000', '10000', '30000', ''), # '' represents all regions
    NNScore = 0
  ) %>%
    # filter(!(Method == 'smwd' & NTopRegions == '30000')) %>%
    filter(!(Method == 'scmet' & NTopRegions %in% c('10000', '30000'))) %>%
    # filter(!(Method == 'scmet' & NTopRegions %in% c('10000'))) %>%
    add_row(Method = 'vseq_cr', NTopRegions = '')
  for (i in 1:nrow(score.df)){
    score.df$NNScore[i] <- nnScoreMethod(
      method = score.df$Method[i], 
      top_n_regions = score.df$NTopRegions[i]
    )
    print(i)
  }
  score.df <- score.df %>% mutate(NTopRegions = as.integer(NTopRegions))
  
  ## Add in total number of regions to summary (score.df)
  res_region <- list(
    vseq = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs")),
    vseq_cr = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs")),
    scbs = loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs")),
    smwd = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs")),
    scmet = loadHDF5SummarizedExperiment(paste0(read_dir, "scmet_regionSummary_vmrs"))
  )
  total_n_regions <- sapply(res_region, function(se) granges(se) %>% length())
  for (i in which(is.na(score.df$NTopRegions))){
    score.df$NTopRegions[i] = total_n_regions[score.df$Method[i]]
  }
  score.df <- score.df %>% mutate(Method = map_chr(Method, methodName))
  
  fwrite(score.df, here(write_dir, paste0("nearest_neighbor_score_", cell_type, "CellType_k", k, "_theta", theta, ".csv")))
}

# computeScoreForAll(k = 100, theta = 0.7, cell_type = 'broad')
# computeScoreForAll(k = 100, theta = 0.7, cell_type = 'sub')
# computeScoreForAll(k = 100, theta = 0.9, cell_type = 'broad')
computeScoreForAll(k = 100, theta = 0.9, cell_type = 'sub')

# ---- Plotting ----

nnScorePlot <- function(k, theta, ylim, ybreaks) {
  
  score_broad.df <- fread(here(write_dir, paste0("nearest_neighbor_score_broadCellType_k", k, "_theta", theta, ".csv")))
  score_sub.df   <- fread(here(write_dir, paste0("nearest_neighbor_score_subCellType_k", k, "_theta", theta, ".csv")))
  score.df <- rbind(data.frame(score_broad.df, Label = 'broad types'),
                    data.frame(score_sub.df,   Label = 'subtypes'))
  p <- score.df %>%
    filter(Method != 'vmrseq CRs') %>%
    ggplot(aes(x = NTopRegions, y = NNScore, color = Method, shape = Label, linetype = Label)) +
    geom_point(size = 2) +
    geom_path() + 
    geom_point(data = score.df %>% filter(Method == 'vmrseq CRs'), 
               aes(x = NTopRegions, y = NNScore, color = Method, shape = Label), 
               size = 2) +
    scale_color_manual(values = COLORVALUES) +
    scale_shape_manual(values = c(1, 16)) +
    scale_linetype_manual(values = c(2, 1)) +
    scale_x_log10(breaks = c(300, 1000, 3000, 10000, 30000, 100000), labels = scales::comma) +
    scale_y_continuous(breaks = ybreaks, limits = ylim) +
    xlab("N Top Regions") + ylab("Nearest Neighbor Score") +
    guides(shape = guide_legend(title = 'Labeled by'),
           linetype = guide_legend(title = 'Labeled by')) +
    theme_classic()
  p
  ggsave(paste0(plot_dir, "point_nnScore_vs_nTopRgions_k",k,"_theta",theta,".png"), height = 5, width = 6)
  ggsave(paste0(plot_dir, "point_nnScore_vs_nTopRgions_k",k,"_theta",theta,"_downsized.png"), height = 3.5, width = 5)
  
  return(p)
}

nnScorePlot(k = 100, theta = 0.7, ylim = c(0, 1), ybreaks = seq(0, 1, 0.5))
nnScorePlot(k = 100, theta = 0.9, ylim = c(0, 1), ybreaks = seq(0, 1, 0.5))


