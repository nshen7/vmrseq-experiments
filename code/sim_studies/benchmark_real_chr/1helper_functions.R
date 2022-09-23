.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
# here::i_am("code/sim_studies/helper_functions/helper_functions_generate_chr.R")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
devtools::load_all("../vmrseq-package/vmrseq/")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(recommenderlab))
## !!! Important: the M_mat stored in a format of recommenderlab::dropNA() 
## outputs, i.e., replace NA as zero and then represent 0 as a very small 
## positive value close to 0.

simPseudoChr <- function(
    N, # total number of cells
    NP, # total number of subpopulations 
    NV, # number of VMRs
    sparseLevel, # sparsity level, should be 1, 2, or 3
    seed, # random seed
    out_dir
) {
  
  set.seed(seed)
  subtype <- "IT-L23_Cux1"
  chromosome <- "chr1"
  
  folder <- paste0("data/interim/sim_studies/real/liu2021_",subtype,"_",chromosome,"/")
  pos0 <- fread(paste0(folder, "pos0.txt.gz")) %>% unlist %>% unname
  # all_file_dirs <- paste0(folder, grep("cell", list.files(folder), value = TRUE))
  
  # Metadata: Proportion of covered CpGs in individual cells
  md <- fread("data/interim/sim_studies/real/liu2021_IT-L23_Cux1_chr1_metadata.csv")
  
  # Sample from designated sparse level
  qt <- quantile(md$cpgCovProp, prob = seq(0.02,0.98,0.32))
  filtered_file_dirs <- with(
    md, 
    fileDir[which(cpgCovProp >= qt[sparseLevel] & cpgCovProp < qt[sparseLevel+1])]
  )
  file_dirs <- sample(filtered_file_dirs, N)
  
  # Read in individual-cell info
  M_mat0 <- NULL
  for (i in 1:N) {
    M_mat0 <- cbind(
      M_mat0,
      fread(file_dirs[i]) %>% as.matrix() %>% dropNA()
    )
    cat(i, " ")
  }
  # M_mat0 <- matrix(NA, nrow = length(pos0), ncol = N)
  # for (i in 1:N) {
  #   M_mat0[, i] <- fread(file_dirs[i]) %>% unlist %>% unname
  #   cat(i, " ")
  # }

  gr <- GRanges(seqnames = "pseudo", ranges = IRanges(start = pos0, end = pos0))
  
  meth <- rowSums(M_mat0)
  total <- rowSums(M_mat0 > 0)
  # meth <- rowSums(M_mat0, na.rm = TRUE)
  # total <- rowSums(!is.na(M_mat0))
  mf <- meth / total

  # Remove sites with 0 coverage 
  index_rm <- which(total <= 0)
  gr <- gr[-index_rm]
  meth <- meth[-index_rm]
  total <- total[-index_rm]
  mf <- mf[-index_rm]
  M_mat0 <- M_mat0[-index_rm, ]
  if (length(gr) != nrow(M_mat0)) stop()
  message("Finished sampling cells.")
  
  # Get CpG clusters from all available CpGs in the subtype pseudo-bulk data (for VMR sampling)
  cluster <- bumphunter::boundedClusterMaker(chr = rep(chromosome, length(pos0)), 
                                             pos = pos0, 
                                             maxGap = 500, 
                                             maxClusterWidth = 10000)
  cluster <- cluster[-index_rm]
  Indexes <- split(seq(along = cluster), cluster)
  lns <- lengths(Indexes)
  Indexes <- Indexes[lns >= 5 & lns <= 500]
  if (length(Indexes) < NV) stop("'NV' is too large.")
  
  # Sample regions with intermediate methylation values preferentially
  mf_means <- map_dbl(Indexes, ~ mean(mf[.x], na.rm = TRUE)) %>% unname()
  # mf_means <- map_dbl(Indexes, ~ mean(colMeans(M_mat0[.x,], na.rm = TRUE), na.rm = TRUE)) %>% unname()
  vmrs_i <- sample(seq_len(length(Indexes)), NV, replace = FALSE, 
                   prob = pmax(1 - sqrt(2) * abs(0.5 - mf_means)^0.5, 0)) %>% sort
  vmrs_Ind <- Indexes[vmrs_i]
  
  # Merge adjacent VMRs if there are equal or less than `maxNumMerge` CpGs in between
  maxNumMerge <- 2
  if (maxNumMerge > 0) {
    for (i in 1:(length(vmrs_Ind)-1)) {
      fr <- vmrs_Ind[[i]]
      bh <- vmrs_Ind[[i+1]]
      if (bh[1] - fr[length(fr)] <= maxNumMerge + 1) {
        combined <- (fr[1]):(bh[length(bh)])
        vmrs_Ind[[i]] <- NA
        vmrs_Ind[[i+1]] <- combined
      }
    }
    vmrs_Ind <- vmrs_Ind[!is.na(vmrs_Ind)]
  }
  print(length(vmrs_Ind)); print(quantile(lengths(vmrs_Ind)))
  message("Finished sampling CpG clusters.")
  
  # Infer indices of null segments in between VMRs
  null_Ind <- list()
  null_Ind[[1]] <- 1:(vmrs_Ind[[1]][1]-1); 
  for (i in 1:length(vmrs_Ind)) {
    h <- vmrs_Ind[[i]][length(vmrs_Ind[[i]])]+1
    r <- ifelse(i < length(vmrs_Ind), yes = vmrs_Ind[[i+1]][1]-1, no = length(total))
    if (h <= r) null_Ind[[i+1]] <- h:r 
  }
  
  # Sample prevalences for VMRs
  pi_prob <- dbinom(1:(NP-1), NP, p = 0.5)
  pi_space <- (1:(NP-1))*1/NP
  pi1s <- sample(pi_space, size = NV, replace = TRUE, prob = pi_prob)
  message("Finished sampling prevalences.")
  
  # Obtain beta prior parameters
  par_u <- vmrseq:::.priorParams(med_cov = median(total), type = "u")
  par_m <- vmrseq:::.priorParams(med_cov = median(total), type = "m")
  pars <- list("u" = par_u, "m" = par_m)
  
  # Fill in single-cell methyl values for VMRs
  M_mat <- NULL
  if (null_Ind[[1]][1] == 1) M_mat <- rbind(M_mat, M_mat0[null_Ind[[1]], ])
  for (i in 1:length(vmrs_Ind)) {
    # rbind vmr
    ix <- vmrs_Ind[[i]]
    miss_mat <- matrix(data = NA, nrow = length(ix), ncol = N)
    miss_mat[as.matrix(M_mat0[ix, ]==0)] <- 1
    vmr_mat <- .sampScMeth2Grp(pi1 = pi1s[i], miss_mat, 0, pars) %>% dropNA()
    M_mat <- rbind(M_mat, vmr_mat)
    # rbind null segment
    if (i+1 <= length(null_Ind)) M_mat <- rbind(M_mat, M_mat0[null_Ind[[i+1]], ])
    cat(i, " ")
  }
  message("Finished filling in single-cell methyl values for VMRs.")
  
  # Useful metadata for sites
  gr$is_vml <- rep(FALSE, length(gr))
  gr$vmr_name <- rep(NA, length(gr))
  gr$pi1 <- rep(NA, length(gr))
  for (i in 1:length(vmrs_Ind)) {
    ix <- vmrs_Ind[[i]]
    gr$is_vml[ix] <- TRUE
    gr$vmr_name[ix] <- i
    gr$pi1[ix] <- pi1s[i]
  }
  gr$meth <- rowSums(M_mat)
  gr$total <- rowSums(M_mat > 0)
  # gr$meth <- rowSums(M_mat, na.rm = TRUE)
  # gr$total <- rowSums(!is.na(M_mat))
  gr$mf <-   gr$meth / gr$total
  message("Finished adding metadata for CpGs.")
  
  write_dir <- paste0(out_dir, "pseudoChr_",
                      subtype, "_", chromosome, "_", 
                      N, "cells_", NP, "subpops_", 
                      NV, "VMRs_sparseLevel", sparseLevel, 
                      "_seed", seed)
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




