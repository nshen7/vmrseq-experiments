source("code/SETPATHS.R")
# devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/cell_type_regional_methylation_scbs/"
if (!file.exists(plot_dir)) dir.create(plot_dir)
md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")
se <- loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs"))
colData(se) <- DataFrame(md)

plotCellTypeMeth <- function(top_i) {
  # index of the top ith region in se
  idx <- order(granges(se)$mcols.var, decreasing = T)[top_i]
  gr <- granges(se)[idx]
  # summarize phenotypic info and methylation value per cell type
  df <- data.frame(
    colData(se),
    MF = assays(se)$M[idx,] / assays(se)$Cov[idx,]
  ) %>% 
    filter(!is.na(MF))
  
  main_title <- paste0("No.", top_i," Top Region")
  sub_title <- paste0(seqnames(gr), ": ", 
                      format(start(gr), big.mark=","), "-", 
                      format(end(gr), big.mark=","),
                      " (width = ", format(end(gr)-start(gr), big.mark=","), ")")
  df %>% ggplot(aes(Neuron.type, MF, fill = Neuron_type1)) +
    geom_boxplot(color = "darkgrey", alpha = 0.8) +
    geom_jitter(color = "darkgrey", height = 0.02, alpha = 0.3, size = 0.5) +
    scale_fill_brewer(palette="Set2") + 
    labs(x = "Cell Type Annotated by Luo at al. (2017)",
         y = "Regional Average Methylation",
         fill = "Neuron Type") +
    ggtitle(main_title, subtitle = sub_title) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave(paste0(plot_dir, "region_topNo", top_i, ".png"), width = 6, height = 3.5)
} 

for (i in 1:100){
  plotCellTypeMeth(i)
  cat(i, " ")
}
