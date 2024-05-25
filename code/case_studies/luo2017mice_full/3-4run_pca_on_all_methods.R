# Ref: SINBAD: a flexible tool for single cell DNA methylation data
#      Yasin Uzun, Wenbao Yu, Changya Chen, Kai Tan
#      bioRxiv 2021.10.23.465577; doi: https://doi.org/10.1101/2021.10.23.465577

source("code/SETPATHS.R")
read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary"
write_dir <- "data/interim/case_studies/luo2017mice_full/pca_on_all_methods"
if (!file.exists(write_dir)) dir.create(write_dir)

# ---- utils -----

imputeAndRunPCA <- function(mf_mat, n_pcs) {
  # `mf_mat`: dimension = n_regions x n_cells
  
  # Imputation suggested by SINBAD paper: 
  # For each feature, use the population mean to replace missing values;
  # i.e., for each row, use the row mean to impute.
  imputed_mf_mat <- mf_mat
  for(i in 1:nrow(imputed_mf_mat)) {
    x <- imputed_mf_mat[i, ]
    imputed_mf_mat[i, is.na(x)] <- mean(x, na.rm = TRUE)
  }
  
  # Run PCA
  pca_result <- prcomp(imputed_mf_mat, scale = TRUE, center = TRUE)
  
  return(pca_result$rotation[, 1:n_pcs])
}

# ---- impute and PCA on regions from each method ----

res_region <- list(
  vseq = loadHDF5SummarizedExperiment(here(read_dir, "vmrseq_regionSummary_vmrs")),
  vseq_cr = loadHDF5SummarizedExperiment(here(read_dir, "vmrseq_regionSummary_crs")),
  scbs = loadHDF5SummarizedExperiment(here(read_dir, "scbs_regionSummary_vmrs")),
  smwd = loadHDF5SummarizedExperiment(here(read_dir, "smallwood_regionSummary_vmrs")),
  scmet = loadHDF5SummarizedExperiment(here(read_dir, "scmet_regionSummary_vmrs"))
)

wrapper <- function(method, n_pcs, top_n_regions = NULL) {
  
  metric <- switch(method,
                   "vseq" = granges(res_region[[method]])$loglik_diff,
                   "scbs" = granges(res_region[[method]])$mcols.var,
                   "smwd" = granges(res_region[[method]])$var_lb,
                   "scmet" = granges(res_region[[method]])$tail_prob)
  if (!is.null(top_n_regions) & method=="vseq_cr") stop("CRs from vmrseq does not have rank.")
  if (!is.null(top_n_regions) & method=="100kbins") stop("100kb bins does not have rank.")
  
  if (!is.null(top_n_regions)) top_ind <- order(metric, decreasing = TRUE)[1:top_n_regions]
  name_seg <- ifelse(is.null(top_n_regions), yes = "", no = paste0("_top", top_n_regions, "regions"))
  
  if (method == '100kbins') {
    
    # Read in matrix of fractional methylation 
    mf_mat0 <- fread(here(read_dir, "..", "100kbins", "input", "100kbins_input.txt.gz")) %>% as.matrix() # dim = 24022 3069
    
    # # Prepocessing: remove features/rows that has low cell coverage
    # n_covered_cells <- apply(mf_mat0, 1, function(x) sum(!is.na(x)))
    # mf_mat <- mf_mat0[n_covered_cells >= 5, ] # dim = 24005  3069
    
    # bins.gr <- do.call(c, readRDS(here(read_dir, "..", "100kbins", "input", "100kbins_feature_metadata.rds")))
    # bins_filtered.gr <- bins.gr[n_covered_cells >= 5]
    # saveRDS(bins_filtered.gr, here(read_dir, "..", "100kbins", "input", "100kbins_feature_metadata_filtered.rds"))
    
  } else if (method %in% c('vseq', 'scbs', 'smwd', 'scmet')) {
    
    se <- res_region[[method]]
    if (!is.null(top_n_regions)) se <- se[top_ind]
    mf_mat <- as.matrix(assays(se)$M/assays(se)$Cov) 
    
  } else if (method == 'vseq_cr') {
    
    se <- res_region[[method]]
    mf_mat <- as.matrix(assays(se)$M/assays(se)$Cov) 
    
  } else {
    stop('Wrong method name!')
  }
  
  # Impute and run PCA
  pca_res <- imputeAndRunPCA(mf_mat, n_pcs = n_pcs)
  fwrite(pca_res, file = here(write_dir, paste0("loadings_", n_pcs, "pcs_", method, name_seg, ".txt.gz")))
}

wrapper(method = '100kbins', n_pcs = 10)
wrapper(method = 'vseq', n_pcs = 10)
wrapper(method = 'vseq_cr', n_pcs = 10)
wrapper(method = 'scbs', n_pcs = 10)
wrapper(method = 'smwd', n_pcs = 10)
wrapper(method = 'scmet', n_pcs = 10)

wrapper(method = 'vseq', n_pcs = 10, top_n_regions = 300)
wrapper(method = 'scbs', n_pcs = 10, top_n_regions = 300)
wrapper(method = 'smwd', n_pcs = 10, top_n_regions = 300)
wrapper(method = 'scmet', n_pcs = 10, top_n_regions = 300)

wrapper(method = 'vseq', n_pcs = 10, top_n_regions = 1000)
wrapper(method = 'scbs', n_pcs = 10, top_n_regions = 1000)
wrapper(method = 'smwd', n_pcs = 10, top_n_regions = 1000)
wrapper(method = 'scmet', n_pcs = 10, top_n_regions = 1000)

wrapper(method = 'vseq', n_pcs = 10, top_n_regions = 3000)
wrapper(method = 'scbs', n_pcs = 10, top_n_regions = 3000)
wrapper(method = 'smwd', n_pcs = 10, top_n_regions = 3000)
wrapper(method = 'scmet', n_pcs = 10, top_n_regions = 3000)

wrapper(method = 'vseq', n_pcs = 10, top_n_regions = 10000)
wrapper(method = 'scbs', n_pcs = 10, top_n_regions = 10000)
wrapper(method = 'smwd', n_pcs = 10, top_n_regions = 10000)

wrapper(method = 'vseq', n_pcs = 10, top_n_regions = 30000)
wrapper(method = 'scbs', n_pcs = 10, top_n_regions = 30000)
wrapper(method = 'smwd', n_pcs = 10, top_n_regions = 30000)
