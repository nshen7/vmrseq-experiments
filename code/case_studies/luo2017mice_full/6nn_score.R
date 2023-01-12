.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
source("code/case_studies/util_functions/6nnScore.R")

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)
md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")

## Parameters for nearest neighbor score
k <- 50; theta <- 0.9

## Nearest neighbor score computed from regional mean methyl for each method
nnScoreMethod <- function(method, top_n_regions = NULL) {
  
  name_seg <- ifelse(top_n_regions == '', yes = '', no = paste0("_top", top_n_regions, "regions"))
  path <- paste0(write_dir, "dissimilarity_matrix_regional_methyl_", method, name_seg, ".txt.gz")
  if (file.exists(path)) d_mat <- fread(path, drop = 1)
  
  stopifnot(all(md$sample == names(d_mat)))
  true_clust <- md$Neuron.type
  return(nnScore(d_mat, true_clust, k, theta))
}

## Summarize nn score from various top n regions
score.df <- expand.grid(
  Method = c('vseq', 'scbs', 'smwd'), 
  NTopRegions = c('300', '1000', '3000', '10000', '20000', ''), # '' represents all regions
  NNScore = 0
) %>%
  filter(!(Method == 'smwd' & NTopRegions == '20000')) %>%
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
  smwd = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs"))#,
  # scmet = loadHDF5SummarizedExperiment(paste0(read_dir, "scmet_regionSummary_vmrs"))
)
total_n_regions <- sapply(res_region, function(se) granges(se) %>% length())
for (i in which(is.na(score.df$NTopRegions))){
  score.df$NTopRegions[i] = total_n_regions[score.df$Method[i]]
}


## Plots
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vseq" = COLORS[1], "vseq_cr" = COLORS[2],
                 "scbs" = COLORS[3], "smwd" = COLORS[4])
# COLORVALUES <- c("vseq" = COLORS[1], "vseq_cr" = COLORS[2],
#                  "scbs" = COLORS[3], "smwd" = COLORS[4], "scmet" = COLORS[5])

score.df %>%
  ggplot(aes(x = NTopRegions, y = NNScore, color = Method)) +
  geom_point() +
  geom_path() + 
  geom_point(data = score.df %>% filter(Method == 'vseq_cr'), 
             aes(x = NTopRegions, y = NNScore, color = Method), 
             size = 3) +
  scale_color_manual(values = COLORVALUES) +
  scale_x_log10() +
  # scale_fill_manual(values = COLORVALUES) +
  xlab("N Top Regions") + ylab("Nearest Neighbor Score") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_line(color = "light grey"))
ggsave(paste0(plot_dir, "point_nnScore_vs_nTopRgions_k",k,"_theta",theta,".png"), height = 4, width = 5)


