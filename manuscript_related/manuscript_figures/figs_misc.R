source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)

plot_dir <- "manuscript_related/manuscript_figures/misc"
if (!file.exists(plot_dir)) dir.create(plot_dir)

# ---- Illustration of between-CpG correlation ----

cell.df <- fread('data/raw_counts/raw_counts_Luo2017_human/Pool_1_AD008_indexed_R1_bismark.tsv.gz') %>%
  mutate(mf = mc_count / total) %>%
  group_by(chr) %>%
  mutate(mf_lag_1 = lag(mf),
         dist_to_lag_1 = c(NA, diff(pos)))

cor.df <- cell.df %>%
  group_by(dist_to_lag_1) %>%
  summarise(cor_mf = cor(mf, mf_lag_1))

cor.df %>%
  ggplot(aes(dist_to_lag_1, cor_mf)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point(size = 0.5, color = 'pink3') +
  xlim(0, 5000) + xlab('Distance of adjacent CpG pair') +
  ylim(-0.5, 1) + ylab('Autocorrelation of adjacent CpG pair') +
  theme_classic()
ggsave(here(plot_dir, 'point_cpg_autocorr_vs_distance_1cellFromLuo.png'), width = 7, height = 5)




