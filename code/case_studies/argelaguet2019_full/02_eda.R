source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(HDF5Array)
library(SummarizedExperiment)

read_dir_rna <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', 'rna')
read_dir_met <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', 'met', 'vmrseq', 'input')
# write_dir <- 
# if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- here('plots', 'case_studies', 'argelaguet2019_full', '02_eda')
if (!file.exists(plot_dir)) dir.create(plot_dir)

md <- fread(here('data', 'metadata', 'argelaguet2019', 'argelaguet2019_full_met&rna_sample_metadata_processed.csv'))

# ---- Misc things ----
## How scMET paper selected the subset of cells
select_idx <- with(md, which(pass_metQC & pass_rnaQC & lineage10x_2 != 'Visceral_endoderm')) 
length(select_idx) # 848

# ---- EDA on cell cycle ----
md %>%
  ggplot(aes(lineage10x_2, G2M.Score)) +
  geom_boxplot() + 
  geom_jitter(color = 'pink', alpha = 0.5, size = 0.7) +
  facet_grid(. ~ stage, scales = "free_x", space = "free_x")
ggsave(here(plot_dir, 'boxplot_cellCycleScore_G2M_vs_lineage10x_2.png'), width = 6, height = 4)

md %>%
  ggplot(aes(lineage10x_2, S.Score)) +
  geom_boxplot() + 
  geom_jitter(color = 'pink', alpha = 0.5, size = 0.7) +
  facet_grid(. ~ stage, scales = "free_x", space = "free_x")
ggsave(here(plot_dir, 'boxplot_cellCycleScore_S_vs_lineage10x_2.png'), width = 6, height = 4)


# ---- EDA on methylation data ----
chr1.se <- loadHDF5SummarizedExperiment(here(read_dir_met, '1'))
dim(chr1.se)
assays(chr1.se)
M_mat <- assays(chr1.se)$M_mat %>% vmrseq::HDF5NAdrop2matrix()

## CpG percentage coverage per cell
cell_covered_cpg <- colSums(M_mat >= 0, na.rm = T) / nrow(M_mat) 
quantile(cell_covered_cpg)
hist(cell_covered_cpg, breaks = 20)

## Global average methylation level
cell_mean_met <- colMeans(M_mat, na.rm = T)

colData(chr1.se) %>%
  as.data.frame() %>%
  mutate(cell_mean_met = cell_mean_met) %>%
  ggplot(aes(lineage10x, cell_mean_met)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(color = 'pink', alpha = 0.5, size = 0.7) +
  facet_grid(. ~ stage, scales = "free_x", space = "free_x") + 
  ylab('Global methylation level per cell') +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
ggsave(here(plot_dir, 'boxplot_cellMeanMet_vs_lineage10x_chr1.png'), width = 6, height = 4)

colData(chr1.se) %>%
  as.data.frame() %>%
  mutate(cell_mean_met = cell_mean_met) %>%
  ggplot(aes(lineage10x_2, cell_mean_met)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(color = 'pink', alpha = 0.5, size = 0.7) +
  facet_grid(. ~ stage, scales = "free_x", space = "free_x") + 
  ylab('Global methylation level per cell') +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) 
ggsave(here(plot_dir, 'boxplot_cellMeanMet_vs_lineage10x_2_chr1.png'), width = 6, height = 4)


