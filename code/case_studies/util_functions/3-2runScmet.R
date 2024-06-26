library(scMET)
library(tidyverse)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)

### Run all chromosomes at one time
runScmet <- function(read_dir, write_dir, plot_dir, efdr, n_cores, min_read = 3, min_cell = 5) {

  dt <- fread(read_dir, col.names = c("Feature", "Cell", "total_reads", "met_reads")) %>%
    mutate(met_reads = met_reads) %>%
    dplyr::filter(total_reads >= min_read) %>%  # exclude regions that have less than 3 CpGs covered
    group_by(Feature) %>%
    dplyr::filter(n() >= min_cell) # exclude features that did not have CpG coverage in at least 5 cells

  # Run model
  fit_obj <- scmet(Y = dt, L = 4, iter = 50000, seed = 12, n_cores = n_cores)

  # plot mean-overdispersion
  gg1 <- scmet_plot_mean_var(obj = fit_obj, y = "gamma",
                             task = NULL, show_fit = TRUE)
  gg2 <- scmet_plot_mean_var(obj = fit_obj, y = "epsilon",
                             task = NULL, show_fit = TRUE)
  cowplot::plot_grid(gg1, gg2, ncol = 2)
  ggsave(paste0(plot_dir, "mean_overdisp_relationship.png"), width = 10, height = 5)

  # store results from various efdr level
  fit_obj <- scmet_hvf(scmet_obj = fit_obj, delta_e = 1-efdr, efdr = efdr)

  fwrite(
    fit_obj$hvf$summary,
    paste0(write_dir, "scmet_summary_output_efdr",efdr,".csv"),
    row.names = FALSE, col.names = TRUE
  )
}


### Run 1 chromosomes at one time 
runScmet1Chr <- function(read_dir, write_dir, plot_dir, chr, n_cores, min_read = 3, min_cell = 5) {
  
  dt <- fread(read_dir, col.names = c("Feature", "Cell", "total_reads", "met_reads")) %>% 
    dplyr::filter(grepl(pattern = paste0(chr,"_"), x = Feature)) %>% # subset to 1 chr
    mutate(met_reads = met_reads) %>%
    dplyr::filter(total_reads >= min_read) %>%  # exclude regions that have less than 3 CpGs covered
    group_by(Feature) %>%
    dplyr::filter(n() >= min_cell) # exclude features that did not have CpG coverage in at least 5 cells
  
  # Run model
  fit_obj <- scmet(Y = dt, L = 4, iter = 50000, seed = 12, n_cores = n_cores)
  saveRDS(fit_obj, paste0(write_dir, "fit_obj_", chr, ".rds"))
  
  # plot mean-overdispersion
  gg1 <- scmet_plot_mean_var(obj = fit_obj, y = "gamma", 
                             task = NULL, show_fit = TRUE)
  gg2 <- scmet_plot_mean_var(obj = fit_obj, y = "epsilon", 
                             task = NULL, show_fit = TRUE)
  cowplot::plot_grid(gg1, gg2, ncol = 2)
  ggsave(paste0(plot_dir, "mean_overdisp_relationship_", chr, ".png"), width = 10, height = 5)
}

hvfScmet <- function(read_dir, write_dir) { 
  
  # combine feature specific estimates for HVF selection
  chrs <- paste0('chr', 1:19)
  hvf_all_chr <- map_dfr(chrs, function(chr) {
    fit_obj <- readRDS(paste0(read_dir, "fit_obj_", chr, ".rds"))
    hvf <- scmet_hvf(scmet_obj = fit_obj)$hvf$summary %>% filter(is_variable) # default settings in scmet
    return(hvf)
  })
  # hvf_all_chr <- hvf_all_chr %>% arrange(desc(gamma))
  
  fwrite(
    hvf_all_chr, 
    paste0(write_dir, "scmet_summary_output_defaultEfdr.csv"), 
    row.names = FALSE, col.names = TRUE
  )
}