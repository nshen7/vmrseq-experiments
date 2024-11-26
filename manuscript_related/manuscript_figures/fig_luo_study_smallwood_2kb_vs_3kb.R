source("code/SETPATHS.R")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
plot_dir <- "manuscript_related/manuscript_figures/fig_luo_study_smallwood_2kb_vs_3kb"
if (!file.exists(plot_dir)) dir.create(plot_dir)

md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")
rownames(md) <- md$sample

res_region <- list(
  smwd = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs")),
  smwd_2kb = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_2kb_regionSummary_vmrs"))
)
methodName <- function(method) switch (method,
                                       'smwd' = 'Smallwood 3kb (default)',
                                       'smwd_2kb' = 'Smallwood 2kb')
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")
COLORVALUES <- c("Smallwood 3kb (default)" = COLORS[5], "Smallwood 2kb" = COLORS[3])
LEVELS <- c('Smallwood 3kb (default)', 'Smallwood 2kb')


# ---- Nearest neighbor scores ----

nnScorePlot <- function(k, theta, ylim, ybreaks) {
  
  score_broad.df <- fread(here(read_dir, paste0("nearest_neighbor_score_broadCellType_k", k, "_theta", theta, ".csv"))) 
  score_sub.df   <- fread(here(read_dir, paste0("nearest_neighbor_score_subCellType_k", k, "_theta", theta, ".csv")))
  score.df <- rbind(data.frame(score_broad.df, Label = 'Broad\nClasses'),
                    data.frame(score_sub.df,   Label = 'Subtypes')) |>
    dplyr::filter(Method %in% c('Smallwood', 'Smallwood 2kb')) |>
    mutate(Method = fct_recode(Method, "Smallwood 3kb (default)" = "Smallwood"))
  score.df %>%
    filter(Method != 'vmrseq CRs') %>%
    ggplot(aes(x = NTopRegions, y = NNScore, color = Method, shape = Label, linetype = Label)) +
    geom_point(size = 3) +
    geom_path() + 
    scale_color_manual(values = COLORVALUES) +
    scale_shape_manual(values = c(1, 16)) +
    scale_linetype_manual(values = c(2, 1)) +
    # scale_x_continuous(labels = scales::comma) +
    scale_x_log10(breaks = c(300, 1000, 3000, 10000, 30000, 100000), labels = scales::comma) +
    scale_y_continuous(breaks = ybreaks, limits = ylim) +
    xlab("# Top Regions") + ylab("Nearest Neighbor Count Score") +
    guides(shape = guide_legend(title = 'Labeled by'),
           linetype = guide_legend(title = 'Labeled by')) +
    theme_classic()
  ggsave(here(plot_dir, paste0("point_nnScore_vs_nTopRgions_k",k,"_theta",theta,".png")), height = 3.5, width = 4.5)
  
  p <- score.df %>%
    filter(Method != 'vmrseq CRs') %>%
    ggplot(aes(x = NSites, y = NNScore, color = Method, shape = Label, linetype = Label)) +
    geom_point(aes(size = NTopRegions)) +
    # geom_point(size = 3) +
    geom_path() +
    scale_size_continuous(
      name = "# top regions",
      breaks = c(300, 1000, 3000, 10000, 30000, 100000)
    ) +
    scale_color_manual(values = COLORVALUES) +
    scale_shape_manual(values = c(1, 16)) +
    scale_linetype_manual(values = c(2, 1)) +
    # scale_x_continuous(labels = scales::comma) +
    scale_x_log10(labels = scales::comma) +
    scale_y_continuous(breaks = ybreaks, limits = ylim) +
    xlab("# CpGs") + ylab("Nearest Neighbor Count Score") +
    guides(shape = guide_legend(title = 'Labeled by'),
           linetype = guide_legend(title = 'Labeled by')) +
    theme_classic()
  p + theme(legend.position = "none")
  ggsave(here(plot_dir, paste0("point_nnScore_vs_nSites_k",k,"_theta",theta,".png")), height = 3.5, width = 3.5)
  
  png(here(plot_dir, "point_nnScore_vs_nSites_subtype_legend.png"), height = 1000, width = 400, res = 200)
  legend <- cowplot::get_legend(p)
  grid::grid.newpage()
  grid::grid.draw(legend)
  dev.off()
}

nnScorePlot(k = 100, theta = 0.7, ylim = c(0, 1), ybreaks = seq(0, 1, 0.5))

# ---- Silhouette scores ----

silhouetteScorePlot <- function() {
  
  score_broad.df <- fread(here(read_dir, paste0("silhouette_score_broadCellType.csv"))) %>% 
    mutate(Label = 'broad types') 
  score_sub.df   <- fread(here(read_dir, paste0("silhouette_score_subCellType.csv"))) %>% 
    mutate(Label = 'subtypes')
  
  score.df <- rbind(data.frame(score_broad.df, Label = 'broad types'),
                    data.frame(score_sub.df,   Label = 'subtypes')) |>
    dplyr::filter(Method %in% c('Smallwood', 'Smallwood 2kb')) |>
    mutate(Method = fct_recode(Method, "Smallwood 3kb (default)" = "Smallwood"))
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
    scale_x_log10(labels = scales::comma) +
    ylim(-0.0551, 0.122) +
    xlab("# Top Regions") + ylab("Silhouette Score") +
    guides(shape = guide_legend(title = 'Labeled by'),
           linetype = guide_legend(title = 'Labeled by')) +
    theme_classic() 
  ggsave(here(plot_dir, "point_silhouetteScore_vs_nTopRgions_subtype.png"), height = 3.5, width = 8)
  
  ## Number of CpG sites in top regions as x axis
  p <- score.df %>%
    filter(Method != 'vmrseq CRs') %>%
    ggplot(aes(x = NSites, y = silhouetteScore, color = Method, shape = Label, linetype = Label)) +
    geom_point(aes(size = NTopRegions)) +
    geom_path() + 
    scale_size_continuous(
      name = "# top regions",
      breaks = c(300, 1000, 3000, 10000, 30000, 100000)
    ) +
    facet_wrap(~Label) +
    scale_color_manual(values = COLORVALUES) +
    scale_shape_manual(values = c(1, 16)) +
    scale_linetype_manual(values = c(2, 1)) +
    scale_x_log10(labels = scales::comma) +
    ylim(-0.0551, 0.122) +
    xlab("# CpGs") + ylab("Silhouette Score") +
    guides(shape = guide_legend(title = 'Labeled by'),
           linetype = guide_legend(title = 'Labeled by')) +
    theme_classic() 
  p + theme(legend.position = "none")
  ggsave(here(plot_dir, "point_silhouetteScore_vs_nSites_subtype.png"), height = 3.5, width = 6.5)
  
  png(here(plot_dir, "point_silhouetteScore_vs_nSites_subtype_legend.png"), height = 1000, width = 400, res = 200)
  legend <- cowplot::get_legend(p)
  grid::grid.newpage()
  grid::grid.draw(legend)
  dev.off()
}

silhouetteScorePlot()


