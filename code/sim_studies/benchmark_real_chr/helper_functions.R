.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
# here::i_am("code/sim_studies/helper_functions/helper_functions_generate_chr.R")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
devtools::load_all("../vmrseq-package/vmrseq/")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(SummarizedExperiment))

simPseudoChr <- function(
    N, # total number of cells
    NP, # total number of subpopulations 
    NV, # number of VMRs
    seed = 1, # random seed
    out_dir
) {
  
  set.seed(seed)
  subtype <- "IT-L23_Cux1"
  chromosome <- "chr1"
  
  folder <- paste0("data/interim/sim_studies/real/liu2021_",subtype,"_",chromosome,"/")
  pos0 <- fread(paste0(folder, "pos0.txt.gz")) %>% unlist %>% unname
  all_file_dirs <- paste0(folder, grep("cell", list.files(folder), value = TRUE))
  file_dirs <- sample(all_file_dirs, N)
  
  # Initialize
  M_mat <- matrix(NA, nrow = length(pos0), ncol = N)
  for (i in 1:N) {
    M_mat[, i] <- fread(file_dirs[i]) %>% unlist %>% unname
    cat(i, " ")
  }
  
  gr <- GRanges(seqnames = "pseudo", ranges = IRanges(start = pos0, end = pos0))
  
  meth <- rowSums(M_mat, na.rm = TRUE)
  total <- rowSums(!is.na(M_mat))
  mf <- meth / total

  # Remove sites with coverage less than 3
  index_rm <- which(total < 3)
  gr <- gr[-index_rm]
  meth <- meth[-index_rm]
  total <- total[-index_rm]
  mf <- mf[-index_rm]
  M_mat <- M_mat[-index_rm, ]
  if (length(gr) != nrow(M_mat)) stop()
  message("Finished sampling cells.")
  
  # Get CpG clusters for VMR sampling
  cluster <- bumphunter::boundedClusterMaker(chr = seqnames(gr), 
                                             pos = start(gr), 
                                             maxGap = 500, 
                                             maxClusterWidth = 10000)
  # cluster <- bumphunter::clusterMaker(chr = seqnames(gr),
  #                                     pos = start(gr),
  #                                     maxGap = 500)
  Indexes <- split(seq(along = cluster), cluster)
  lns <- lengths(Indexes)
  Indexes <- Indexes[lns >= 5 & lns <= 500]
  if (length(Indexes) < NV) stop("'NV' is too large.")
  
  # Sample regions with intermediate methylation values preferentially
  mf_meds <- map_dbl(Indexes, ~ median(mf[.x], na.rm = TRUE)) %>% unname()
  vmrs_i <- sample(seq_len(length(Indexes)), NV, replace = FALSE, 
                   prob = pmax(1 - sqrt(2) * abs(0.5 - mf_meds)^0.5, 0)) %>% sort
  vmrs_Ind <- Indexes[vmrs_i]
  message("Finished sampling CpG clusters.")
  
  # # Sample prevalences for VMRs
  # p_meth <- 0.77
  # pi_space <- (1:(NP-1))/NP
  # pi_prob <- dbinom(1:(NP-1), NP, p = p_meth)
  # pi1s <- sample(pi_space, size = NV, replace = TRUE, prob = pi_prob)
  
  # (Original) Sample prevalences for VMRs
  pi_prob <- dbinom(1:(NP-1), NP, p = 0.5)
  pi_space <- (1:(NP-1))*1/NP
  pi1s <- sample(pi_space, size = NV, replace = TRUE, prob = pi_prob)
  
  # Fill in single-cell methyl values for VMRs
  par_u <- vmrseq:::.priorParams(med_cov = median(total), type = "u")
  par_m <- vmrseq:::.priorParams(med_cov = median(total), type = "m")
  pars <- list("u" = par_u, "m" = par_m)
  
  for (i in 1:NV) {
    ix <- vmrs_Ind[[i]]
    miss_mat <- matrix(data = NA, nrow = length(ix), ncol = N) 
    miss_mat[!is.na(M_mat[ix, ])] <- 1
    M_mat[ix, ] <- .sampScMeth2Grp(pi1 = pi1s[i], miss_mat, 0, pars)
  }
  message("Finished filling in single-cell methyl values for VMRs.")
  
  # Useful metadata for sites
  gr$is_vml <- rep(FALSE, length(gr))
  gr$vmr_name <- rep(NA, length(gr))
  gr$pi1 <- rep(NA, length(gr))
  for (i in 1:NV) {
    ix <- vmrs_Ind[[i]]
    gr$is_vml[ix] <- TRUE
    gr$vmr_name[ix] <- i
    gr$pi1[ix] <- pi1s[i]
  }
  gr$meth <- rowSums(M_mat, na.rm = TRUE)
  gr$total <- rowSums(!is.na(M_mat))
  gr$mf <-   gr$meth /   gr$total
  message("Finished adding metadata for CpGs.")
  
  write_dir <- paste0(out_dir, "pseudoChr_",
                      subtype, "_", chromosome, "_", 
                      N, "cells_", NP, "subpops_", 
                      NV, "VMRs_seed", seed)
  saveHDF5SummarizedExperiment(
    x = SummarizedExperiment(
      assays = list(M_mat = M_mat), 
      rowRanges = gr
    ),
    dir = write_dir,
    replace = TRUE
  )
  
  message(paste0("Finished saving pseudo-chromosome in a HDF5SummarizedExperiment object to:\n '",
                 write_dir, "'"))
  return(write_dir)
}



# ==== utils ====
.sampBetaBinom <- function(bd, state, pars) {
  if (!state%in%c(0,1)) stop("'state' should be either 0 or 1")
  if (!state) {
    par_u <- pars[['u']]
    return(rZIBB(1, mu = par_u['mu'], sigma = par_u['sigma'], nu = par_u['nu'], bd = bd))
  } else {
    par_m <- pars[['m']]
    return(rBB(1, mu = par_m['mu'], sigma = par_m['sigma'], bd = bd))
  }
}

# Sample methylation level for individual cells
.sampMethMat <- function(state_seq, miss_mat, N, sigma, pars) {
  # `sigma` is the noise level
  
  K <- length(state_seq)
  stopifnot("Number of CpGs `K` not equal to length of state sequence." = K == length(state_seq))
  
  mf_mat <- do.call(
    rbind, 
    map(1:K, function(.x) {
      x <- miss_mat[.x, ]
      total <- sum(x, na.rm = T)
      meth <- .sampBetaBinom(bd = total, state = state_seq[.x], pars)
      ind_u <- sample(which(!is.na(x)), total - meth)
      x[ind_u] <- 0
      return(x)
    })
  )
  nois_ind_mat <- do.call(
    rbind,
    replicate(K, rbernoulli(N, p = sigma), simplify = F)
  )
  mf_mat[nois_ind_mat] <- 1 - mf_mat[nois_ind_mat]
  return(mf_mat)
}

# Sample methylation level for individual cells in 2-grouping
.sampScMeth2Grp <- function(pi1, miss_mat, sigma, pars) {

  N <- ncol(miss_mat)
  K <- nrow(miss_mat)
  
  state_seq <- rep(1, K) # ps: hidden states in VMR are all (1,0), which is encoded as 1
  state_vec <- do.call(rbind, map(state_seq, vmrseq:::.translateState2Grp))
  
  ind_g1 <- sample(1:N, ceiling(pi1*N)) %>% sort
  ind_g2 <- which(! 1:N %in% ind_g1)
  mf_mat <- matrix(nrow = K, ncol = N)
  mf_mat[,ind_g1] <- .sampMethMat(state_seq = state_vec[,1], miss_mat = miss_mat[,ind_g1], N = length(ind_g1), sigma, pars)
  mf_mat[,ind_g2] <- .sampMethMat(state_seq = state_vec[,2], miss_mat = miss_mat[,ind_g2], N = length(ind_g2), sigma, pars)
  
  return(mf_mat * miss_mat)
}




