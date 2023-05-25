source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(HDF5Array)
library(SummarizedExperiment)

read_dir_rna <- here('data', 'interim', 'case_studies', 'argelaguet2019', 'rna')
read_dir_met <- here('data', 'interim', 'case_studies', 'argelaguet2019', 'met', 'vmrseq', 'input')
# write_dir <- 
# if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- here('plots', 'case_studies', 'argelaguet2019', '02_eda')
if (!file.exists(plot_dir)) dir.create(plot_dir)

  
# ---- Misc things ----
## How scMET paper selected the subset of cells
select_idx <- with(md, which(pass_metQC & pass_rnaQC & lineage10x_2 != 'Visceral_endoderm')) 
length(select_idx) # 848

# ---- EDA on met data ----
chr1.se <- loadHDF5SummarizedExperiment(here(read_dir_met, '1'))
dim(chr1.se)
assays(chr1.se)
# M_mat <- assays(chr1.se)$M_mat %>% as("sparseMatrix") %>% recommenderlab::dropNA2matrix()
M_mat <- assays(chr1.se)$M_mat %>% vmrseq::HDF5NAdrop2matrix()
cell_mean_met <- colMeans(M_mat, na.rm = T)

colData(chr1.se) %>%
  as.data.frame() %>%
  mutate(cell_mean_met = cell_mean_met) %>%
  ggplot(aes(lineage10x_2, cell_mean_met)) +
  geom_boxplot() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0))
ggsave(here(plot_dir, 'boxplot_cellMeanMet_vs_cellType_chr1.png'), width = 6, height = 4)
