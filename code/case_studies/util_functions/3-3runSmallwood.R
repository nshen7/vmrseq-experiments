source("code/SETPATHS.R")
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
# n_cores <- 22

runSmallwood <- function(read_dir, write_dir, chrs, alpha, n_cores) {
  
  register(MulticoreParam(workers = n_cores))
  
  # load input
  if (any(!chrs %in% list.files(read_dir))) stop("Did not find all chromosomes in `read_dir`")
  SE <- do.call(
    rbind, 
    lapply(paste0(read_dir, chrs), loadHDF5SummarizedExperiment)
  ) 
  
  # run model
  r_hat_ij <- (assays(SE)$M + 1) / (assays(SE)$Cov + 2)
  var_r_hat_ij <- sqrt(r_hat_ij*(1-r_hat_ij) / assays(SE)$Cov)
  
  w_ij <- 1 / var_r_hat_ij
  w_i <- rowSums(w_ij)
  r_bar_i <- rowSums(r_hat_ij * w_ij) / w_i
  
  n_i_sq <- w_i / (w_i^2 - rowSums(w_ij^2))
  v_hat_i <- n_i_sq * rowSums(w_ij * (r_hat_ij - r_bar_i)^2)
  
  n_i <- sqrt(n_i_sq)
  chi_sq_i <- map_dbl(1:length(n_i), ~ qchisq(p = 1-alpha/2, df = n_i[.x]))
  
  v_hat_i_lb <- n_i * v_hat_i / chi_sq_i
  values(SE)$var_lb <- v_hat_i_lb
  
  # save model output
  saveRDS(granges(SE), paste0(write_dir, "smallwood_output_varLowerBound.rds"))
}

