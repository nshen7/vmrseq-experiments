.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
library(scMET)
library(data.table)
library(tidyverse)
library(HDF5Array)
library(GenomicRanges)
library(BiocParallel)
library(DelayedArray)
library(DelayedMatrixStats)
library(SummarizedExperiment)

# register(SnowParam(workers = 6, type = "FORK"))

formatForScmet <- function(se, bp_size) {
  
  cuts <- seq(start(se)[1], start(se)[length(se)], bp_size)
  wds.gr <- GRanges(
    seqnames = seqnames(se)[1],
    ranges = IRanges(start = cuts,
                     end = c(cuts[-1]-1, start(se)[length(se)]))
  )
  hits <- findOverlaps(granges(se), wds.gr)
  
  computeFeature <- function(i) { # i th feature/window
    feat.se <- se[queryHits(hits)[subjectHits(hits)==i],]
    total_reads <- colSums(assays(feat.se)$Cov > 0, na.rm = T)
    met_reads <- colSums(round(assays(feat.se)$M/assays(feat.se)$Cov), na.rm = T)
    feat.df <- data.frame(
      Feature = paste0("Window_", i),
      Cell = paste0("Cell_",1:length(total_reads)),
      total_reads = total_reads,
      met_reads = met_reads
    ) %>% dplyr::filter(total_reads > 0)
    return(feat.df)
  }
  
  feats.df <- do.call(
    rbind,
    bplapply(1:length(wds.gr), computeFeature)
  )
  return(feats.df)
}

se <- loadHDF5SummarizedExperiment(
  "data/processed/processed_liu2021/sample5F_181220_GSE135169_subtype_MSN-D1_Plxnc1_123cells/"
) %>% subset(seqnames=="chr1")
Y <- formatForScmet(se, bp_size = 20000)

# Y <- fread("data/interim/sim_studies/benchmark_pseudo_chr/scmet/input/scmet_input_20kbWindow_1000cells_20subpops.txt.gz")
Y <- Y %>% 
  filter(total_reads >= 5) %>%
  group_by(Feature) %>%
  filter(n() > 5)
hist(Y$met_reads/Y$total_reads)

fit_obj <- scmet(Y = Y, L = 4, iter = 20000, seed = 1)
fit_obj <- scmet_hvf(scmet_obj = fit_obj)
head(fit_obj$hvf$summary)
saveRDS(fit_obj, "code/test/test_scmet/fit_obj.rds")

gg1 <- scmet_plot_mean_var(obj = fit_obj, y = "gamma", 
                           task = NULL, show_fit = TRUE)
gg2 <- scmet_plot_mean_var(obj = fit_obj, y = "epsilon", 
                           task = NULL, show_fit = TRUE)
cowplot::plot_grid(gg1, gg2, ncol = 2)
ggsave("code/test/test_scmet/scmet_plot.png", width = 10, height = 5)

gg1 <- scmet_plot_vf_tail_prob(obj = fit_obj, x = "mu", task = "hvf")
gg2 <- scmet_plot_mean_var(obj = fit_obj, y = "gamma", task = "hvf")
cowplot::plot_grid(gg1, gg2, ncol = 2)
ggsave("code/test/test_scmet/scmet_hvf_plot.png", width = 10, height = 5)









