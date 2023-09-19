source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(pheatmap)

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

res.se <- list(
  'vseq' = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs")),
  'vseq_cr' = loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs")),
  'scbs' = loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs")),
  'smwd' = loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs")),
  'scmet' = loadHDF5SummarizedExperiment(paste0(read_dir, "scmet_regionSummary_vmrs"))
)


# Read in SE object of regional summary of VMRs 
md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")

# ---- Summarize cell type average regional methylation ----

cellTypeMeanMethyl <- function(se, metric, save_dir) {
  
  idx <- which(!duplicated(colData(se)$Neuron.type))
  types <- colData(se)$Neuron.type[idx]
  types_transmitter <- colData(se)$Neuron_type1[idx]
  
  # (average across cells in each cell type)
  typeMeanMeth <- function(type) {
    cell_idx <- which(colData(se)$Neuron.type == type)
    M_mat <- as.matrix(assays(se)$M[, cell_idx] / assays(se)$Cov[, cell_idx])
    n_avail_cell <- apply(M_mat, MARGIN = 1, FUN = function(x) sum(!is.na(x)))
    M_mean <- apply(M_mat, MARGIN = 1, FUN = function(x) mean(x, na.rm = T))
    
    M_mat_other <- as.matrix(assays(se)$M[, -cell_idx] / assays(se)$Cov[, -cell_idx])
    M_mean_other <- apply(M_mat_other, MARGIN = 1, FUN = function(x) mean(x, na.rm = T))
    return(list(N = n_avail_cell, M_mean = M_mean, M_mean_other = M_mean_other))
  }
  types_smr <- lapply(X = types, FUN = typeMeanMeth)
  
  # Matrix of number of non-NA cells in each region (n_region x n_celltype)
  types_N <- do.call(cbind, lapply(types_smr, function(x) x$N)) %>% data.frame()
  colnames(types_N) <- types
  # Matrix of mean regional methylation of each cell type (n_region x n_celltype)
  types_mean <- do.call(cbind, lapply(types_smr, function(x) x$M_mean)) %>% data.frame()
  colnames(types_mean) <- types
  
  types_smr.se <- SummarizedExperiment(assays = SimpleList(reginal_methyl = types_mean, n_cell = types_N),
                                       rowRanges = rowRanges(se),
                                       colData = DataFrame(neuron_type = types_transmitter))
  
  saveRDS(types_smr.se, save_dir)
}

# Run
for (method in names(res.se)) {

  se <- res.se[[method]]
  colData(se) <- DataFrame(md)

  if (method != "vseq_cr") {
    metric <- switch(method,
                     "vseq" = granges(se)$loglik_diff,
                     "scbs" = granges(se)$mcols.var,
                     "smwd" = granges(se)$var_lb,
                     "scmet" = granges(se)$gamma)
  }

  save_dir <- here(write_dir, paste0('SummarizedExperiment_cellType_regionalMean_', method, '.rds'))
  cellTypeMeanMethyl(se = se, metric = metric, save_dir = save_dir)
}

# ---- Shuffle cell type labels and summarize fake cell type average regional methylation ----


shuffledCellTypeMeanMethyl <- function(se, metric, save_dir, seed) {
  
  set.seed(seed = seed)
  shuffled_type <- sample(colData(se)$Neuron.type)
  
  idx <- which(!duplicated(shuffled_type))
  types <- shuffled_type[idx]

  # (average across cells in each cell type)
  typeMeanMeth <- function(type) {
    cell_idx <- which(shuffled_type == type)
    M_mat <- as.matrix(assays(se)$M[, cell_idx] / assays(se)$Cov[, cell_idx])
    n_avail_cell <- apply(M_mat, MARGIN = 1, FUN = function(x) sum(!is.na(x)))
    M_mean <- apply(M_mat, MARGIN = 1, FUN = function(x) mean(x, na.rm = T))
    
    # M_mat_other <- as.matrix(assays(se)$M[, -cell_idx] / assays(se)$Cov[, -cell_idx])
    # M_mean_other <- apply(M_mat_other, MARGIN = 1, FUN = function(x) mean(x, na.rm = T))
    # return(list(N = n_avail_cell, M_mean = M_mean, M_mean_other = M_mean_other))
    return(list(N = n_avail_cell, M_mean = M_mean))
  }
  types_smr <- lapply(X = types, FUN = typeMeanMeth)
  
  # Matrix of number of non-NA cells in each region (n_region x n_celltype)
  types_N <- do.call(cbind, lapply(types_smr, function(x) x$N)) %>% data.frame()
  colnames(types_N) <- types
  # Matrix of mean regional methylation of each cell type (n_region x n_celltype)
  types_mean <- do.call(cbind, lapply(types_smr, function(x) x$M_mean)) %>% data.frame()
  colnames(types_mean) <- types
  
  types_smr.se <- SummarizedExperiment(assays = SimpleList(reginal_methyl = types_mean, n_cell = types_N),
                                       rowRanges = rowRanges(se))
  saveRDS(types_smr.se, save_dir)
}


## Run
write_dir_2 <- here(write_dir, 'cell_type_shuffled')
if (!file.exists(write_dir_2)) dir.create(write_dir_2)

for (method in names(res.se)) {
  se <- res.se[[method]]
  colData(se) <- DataFrame(md)

  if (method != "vseq_cr") {
    metric <- switch(method,
                     "vseq" = granges(se)$loglik_diff,
                     "scbs" = granges(se)$mcols.var,
                     "smwd" = granges(se)$var_lb,
                     "scmet" = granges(se)$gamma)
  }

  for (seed in 1:10) {
    save_dir <- here(write_dir_2, paste0('SummarizedExperiment_shuffledCellType_regionalMean_', method, '_seed', seed, '.rds'))
    shuffledCellTypeMeanMethyl(se = se, metric = metric, save_dir = save_dir, seed = seed)
  }
}




