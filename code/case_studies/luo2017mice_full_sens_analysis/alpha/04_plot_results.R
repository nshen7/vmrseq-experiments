source("code/SETPATHS.R")

read_dir <- "data/interim/case_studies/luo2017mice_full_sens_analysis/result_summary/alpha"
write_dir <- "data/interim/case_studies/luo2017mice_full_sens_analysis/result_summary/alpha"
if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- "plots/case_studies/luo2017mice_full_sens_analysis"
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = T)

# ---- utils ----
# Color settings
methodName <- function(method) switch (method,
                                       'vmrs' = 'vmrseq',
                                       'crs' = 'vmrseq CRs')
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = COLORS[1], "vmrseq CRs" = COLORS[2])
k <- 100
theta <- 0.7


# ---- overlap upset plot ----


res_vmrs <- list(
  `alpha = 0.01` = here(read_dir, "alpha0.01/vmrseq_regionSummary_vmrs") |> loadHDF5SummarizedExperiment() |> granges(),
  `alpha = 0.025` = here(read_dir, "alpha0.025/vmrseq_regionSummary_vmrs") |> loadHDF5SummarizedExperiment() |> granges(),
  `alpha = 0.05` = here(read_dir, "alpha0.05/vmrseq_regionSummary_vmrs") |> loadHDF5SummarizedExperiment() |> granges(),
  `alpha = 0.1` = here(read_dir, "alpha0.1/vmrseq_regionSummary_vmrs") |> loadHDF5SummarizedExperiment() |> granges()
)
res_crs <- list(
  `alpha = 0.01` = here(read_dir, "alpha0.01/vmrseq_regionSummary_crs") |> loadHDF5SummarizedExperiment() |> granges(),
  `alpha = 0.025` = here(read_dir, "alpha0.025/vmrseq_regionSummary_crs") |> loadHDF5SummarizedExperiment() |> granges(),
  `alpha = 0.05` = here(read_dir, "alpha0.05/vmrseq_regionSummary_crs") |> loadHDF5SummarizedExperiment() |> granges(),
  `alpha = 0.1` = here(read_dir, "alpha0.1/vmrseq_regionSummary_crs") |> loadHDF5SummarizedExperiment() |> granges()
)

# Overlap between detected regions
comb_mat <- ComplexHeatmap::make_comb_mat(res_vmrs)
png(here(plot_dir, "UpSetPlot_alpha_VMRs_overlapWindowSize.png"), height = 450, width = 600, res = 150)
ComplexHeatmap::UpSet(comb_mat, set_order = 1:4, comb_order = order(ComplexHeatmap::comb_size(comb_mat), decreasing = T))
dev.off()

comb_mat <- ComplexHeatmap::make_comb_mat(res_crs)
png(here(plot_dir, "UpSetPlot_alpha_CRs_overlapWindowSize.png"), height = 450, width = 600, res = 150)
ComplexHeatmap::UpSet(comb_mat, set_order = 1:4, comb_order = order(ComplexHeatmap::comb_size(comb_mat), decreasing = T))
dev.off()

# ---- plot ----
score.df <- fread(here(write_dir, paste0("nearest_neighbor_score_k", k, "_theta", theta, ".csv")))
score.df %>%
  pivot_longer(cols = -1, 
               names_sep = '_', names_prefix = 'NNScore_', names_to = c('Method', 'Label'), 
               values_to = 'NNScore') %>% 
  mutate(Method = map_chr(Method, methodName),
         Label = factor(Label, labels = c('broad' = 'Broad\nClasses', 'sub' = 'Subtypes'))) %>%
  ggplot(aes(alpha, NNScore, color = Method, shape = Label, linetype = Label)) + 
  geom_point(size = 3) + 
  geom_line() +
  scale_color_manual(values = COLORVALUES) +
  scale_shape_manual(values = c(1, 16)) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_x_continuous(breaks = score.df$alpha) +
  scale_y_continuous(breaks = seq(0, 1, 0.5), lim = c(0, 1)) +
  xlab("Alpha Value") + ylab("Nearest Neighbor Count Score") +
  guides(shape = guide_legend(title = 'Labeled by'),
         linetype = guide_legend(title = 'Labeled by')) +
  theme_classic()
ggsave(here(plot_dir, paste0("point_nnScore_vs_alpha_k",k,"_theta",theta,".png")), height = 3.5, width = 4.5)

