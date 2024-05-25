source("code/SETPATHS.R")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(here)

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)
md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")

# ---- utils ----

res_region <- list(
  vseq = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs")),
  vseq_cr = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs")),
  scbs = loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs")),
  smwd = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs")),
  scmet = loadHDF5SummarizedExperiment(paste0(read_dir, "scmet_regionSummary_vmrs"))
)
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

computeScoreForAll <- function(cell_type) {
  
  
  if (cell_type == 'broad') {
    # Use main neuron types as reference cell labels: 'Excitatory' or 'inhibitory'
    true_clust <- md$Neuron_type1
  } else if (cell_type == 'sub') {
    # Use neuron subtypes as reference cell labels: 'mL5-1' etc
    true_clust <- md$Neuron.type
  }
  
  ## silhouette score computed from regional mean methyl for each method
  silhouetteScoreMethod <- function(method, top_n_regions = NULL) {
    
    name_seg <- ifelse(top_n_regions == '', yes = '', no = paste0("_top", top_n_regions, "regions"))
    path <- paste0(read_dir, "dissimilarity_matrix_regional_methyl_", method, name_seg, ".txt.gz")
    if (file.exists(path)) d_mat <- fread(path, drop = 1) else stop('Dissimilarity matrix not written yet!')
    
    stopifnot(all(md$sample == names(d_mat)))
    sil.df <- cluster::silhouette(x = as.integer(factor(true_clust)), dmatrix = as.matrix(d_mat)) %>% as.data.frame()
    return(mean(sil.df$sil_width))
  }
  
  ## Count the number sites in top regions
  countNSites <- function(method, top_n_regions = NULL) {
    
    metric <- switch(method,
                     "vseq" = granges(res_region[[method]])$loglik_diff,
                     "scbs" = granges(res_region[[method]])$mcols.var,
                     "smwd" = granges(res_region[[method]])$var_lb,
                     "scmet" = granges(res_region[[method]])$tail_prob)
    top_ind <- order(metric, decreasing = TRUE)[1:top_n_regions]
    top.gr <- granges(res_region[[method]])[top_ind]
    
    n_sites <- findOverlaps(GenomicRanges::reduce(top.gr), sites.gr) %>% length()
    return(n_sites)
    
  }
  
  ## Summarize silhouette score from various top n regions
  score.df <- expand.grid(
    Method = c('vseq', 'scbs', 'smwd', 'scmet'), 
    NTopRegions = c('300', '1000', '3000', '10000', '30000', ''), # '' represents all regions
    silhouetteScore = 0
  ) %>%
    # filter(!(Method == 'smwd' & NTopRegions == '30000')) %>%
    filter(!(Method == 'scmet' & NTopRegions %in% c('10000', '30000'))) %>%
    # filter(!(Method == 'scmet' & NTopRegions %in% c('10000'))) %>%
    add_row(Method = 'vseq_cr', NTopRegions = '') 
  for (i in 1:nrow(score.df)){
    ## Compute silhouette score
    score.df$silhouetteScore[i] <- silhouetteScoreMethod(
      method = score.df$Method[i], 
      top_n_regions = score.df$NTopRegions[i]
    )
    
    ## Add in number of sites in detected regions
    if (score.df$NTopRegions[i] != "") 
      score.df$NSites[i] <- countNSites(
        method = score.df$Method[i], 
        top_n_regions = score.df$NTopRegions[i]
      )
    
    print(i)
  }
  score.df <- score.df %>% mutate(NTopRegions = as.integer(NTopRegions))
  
  ## Add in total number of regions&sites to summary (score.df)
  total_n_regions <- sapply(res_region, function(se) granges(se) %>% length())
  total_n_sites <- sapply(res_region, function(se) findOverlaps(GenomicRanges::reduce(granges(se)), sites.gr) %>% length())
  for (i in which(is.na(score.df$NTopRegions))){
    score.df$NTopRegions[i] = total_n_regions[score.df$Method[i]]
    score.df$NSites[i]      = total_n_sites  [score.df$Method[i]]
  }
  score.df <- score.df %>% mutate(Method = map_chr(Method, methodName))
  
  fwrite(score.df, here(write_dir, paste0("silhouette_score_", cell_type, "CellType.csv")))
}

# computeScoreForAll(cell_type = 'broad')
# computeScoreForAll(cell_type = 'sub')

# ---- Plotting ----

silhouetteScorePlot <- function() {
  
  score_broad.df <- fread(here(write_dir, paste0("silhouette_score_broadCellType.csv"))) %>% mutate(Label = 'broad types') 
  score_sub.df   <- fread(here(write_dir, paste0("silhouette_score_subCellType.csv"))) %>% mutate(Label = 'subtypes')
  
  score.df <- rbind(data.frame(score_broad.df, Label = 'broad types'),
                    data.frame(score_sub.df,   Label = 'subtypes'))
  ## Number of top regions as x axis
  score.df %>%
    filter(Method != 'vmrseq CRs') %>%
    ggplot(aes(x = NTopRegions, y = silhouetteScore, color = Method, shape = Label, linetype = Label)) +
    geom_point(size = 3) +
    geom_path() + 
    geom_point(data = score.df %>% filter(Method == 'vmrseq CRs'), 
               aes(x = NTopRegions, y = silhouetteScore, color = Method, shape = Label), 
               size = 3) +
    facet_wrap(~Label) +
    scale_color_manual(values = COLORVALUES) +
    scale_shape_manual(values = c(1, 16)) +
    scale_linetype_manual(values = c(2, 1)) +
    # scale_x_continuous(labels = scales::comma) +
    scale_x_log10(labels = scales::comma) +
    # scale_x_log10(breaks = c(300, 1000, 3000, 10000, 30000, 100000), labels = scales::comma) +
    # scale_y_continuous(breaks = ybreaks, limits = ylim) +
    xlab("# Top Regions") + ylab("Silhouette Score") +
    guides(shape = guide_legend(title = 'Labeled by'),
           linetype = guide_legend(title = 'Labeled by')) +
    theme_classic()
  ggsave(paste0(plot_dir, "point_silhouetteScore_vs_nTopRgions_subtype.png"), height = 3.5, width = 8)
  
  ## Number of CpG sites in top regions as x axis
  p <- score.df %>%
    filter(Method != 'vmrseq CRs') %>%
    ggplot(aes(x = NSites, y = silhouetteScore, color = Method, shape = Label, linetype = Label)) +
    geom_point(aes(size = NTopRegions)) +
    # geom_point(size = 3) +
    geom_path() + 
    geom_point(data = score.df %>% filter(Method == 'vmrseq CRs'), 
               # aes(x = NSites, y = silhouetteScore, color = Method, shape = Label), size = 3) +
               aes(x = NSites, y = silhouetteScore, color = Method, shape = Label, size = NTopRegions)) +
    scale_size_continuous(
      name = "# top regions",
      breaks = c(300, 1000, 3000, 10000, 30000, 100000)
    ) +
    facet_wrap(~Label) +
    scale_color_manual(values = COLORVALUES) +
    scale_shape_manual(values = c(1, 16)) +
    scale_linetype_manual(values = c(2, 1)) +
    # scale_x_continuous(labels = scales::comma) +
    scale_x_log10(labels = scales::comma) +
    # scale_y_continuous(breaks = ybreaks, limits = ylim) +
    xlab("# CpGs") + ylab("Silhouette Score") +
    guides(shape = guide_legend(title = 'Labeled by'),
           linetype = guide_legend(title = 'Labeled by')) +
    theme_classic() 
  p + theme(legend.position = "none")
  ggsave(paste0(plot_dir, "point_silhouetteScore_vs_nSites_subtype.png"), height = 3.5, width = 6.5)
  
  png(paste0(plot_dir, "point_silhouetteScore_vs_nSites_subtype_legend.png"), height = 1000, width = 300, res = 200)
  legend <- cowplot::get_legend(p)
  grid::grid.newpage()
  grid::grid.draw(legend)
  dev.off()
}

silhouetteScorePlot()


