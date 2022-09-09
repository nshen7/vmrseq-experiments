.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
# here::i_am("code/sim_studies/helper_functions/helper_functions_generate_chr.R")
source("code/sim_studies/benchmark_sim_chr/helper_functions_generate_region.R")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(GenomicRanges))

# register(SnowParam(workers = 6, type = "FORK"))

# # Import CpG distances in chr 1 from a real dataset
# dat0 <- fread("data/processed/summarized_liu2021/summarized_subtype_MSN-D1_Khdrbs3_540cells.csv.gz") %>%
#   filter(chr == "chr1") %>%
#   mutate(state = as.numeric(cell_meth/cell_cov >= 0.5)) %>%
#   dplyr::select(c(pos, state))

# ==== Main function for generating  object of a pseudo chromosome ====

simPseudoChr <- function(
    N, # total number of cells
    NV, # number of VMRs
    NP, # total number of subpopulations 
    out_dir, # path to where simulated should be stored
    seed = 2022, # random seed
    sigma = 0 # noise level
) {
  
  set.seed(seed)
  
  subtype <- "IT-L23_Cux1"
  chromosome <- "chr1"
  
  folder <- paste0("data/interim/sim_studies/real/liu2021_",subtype,"_",chromosome,"/")
  pos0 <- fread(paste0(folder, "pos0.txt.gz")) %>% unlist %>% unname
  all_file_dirs <- paste0(folder, grep("cell", list.files(folder), value = TRUE))
  file_dirs <- sample(all_file_dirs, N)
  
  # Get missing distribution from real data
  M_mat <- matrix(NA, nrow = length(pos0), ncol = N)
  for (i in 1:N) {
    M_mat[, i] <- fread(file_dirs[i]) %>% unlist %>% unname
    # cat(i, " ")
  }
  
  gr <- GRanges(seqnames = "pseudo", ranges = IRanges(start = pos0, end = pos0))
  
  meth <- rowSums(M_mat, na.rm = TRUE)
  total <- rowSums(!is.na(M_mat))
  state <- as.numeric(meth / total >= 0.5)
  rm(meth)

  # Remove sites with coverage less than 3
  index_rm <- which(total < 3)
  gr <- gr[-index_rm]
  M_mat <- M_mat[-index_rm, ]
  total <- total[-index_rm]
  state <- state[-index_rm]

  message("Cell/site missing distribution obtained from real data.")
  
  # Get CpG clusters for VMR sampling (minCov >= 3 & 5 <= num_cpg <= 500)
  cluster <- bumphunter::boundedClusterMaker(chr = seqnames(gr), 
                                             pos = start(gr), 
                                             maxGap = 500, 
                                             maxClusterWidth = 10000)
  Indexes <- split(seq(along = cluster), cluster)
  lns <- lengths(Indexes)
  min_covs <- sapply(Indexes, function(.x) min(total[.x])) %>% unname
  Indexes <- Indexes[lns >= 5 & lns <= 500 & min_covs >= 3]
  if (length(Indexes) < NV) stop("'NV' is too large.")
  
  # Sample regions 
  vmrs_i <- sample(seq_len(length(Indexes)), NV, replace = FALSE) %>% sort
  inds_2g <- Indexes[vmrs_i]
  # Merge adjacent VMRs if there are equal or less than `maxNumMerge` CpGs in between
  maxNumMerge <- 2
  if (maxNumMerge > 0) {
    for (i in 1:(length(inds_2g)-1)) {
      fr <- inds_2g[[i]]
      bh <- inds_2g[[i+1]]
      if (bh[1] - fr[length(fr)] <= maxNumMerge + 1) {
        combined <- (fr[1]):(bh[length(bh)])
        inds_2g[[i]] <- NA
        inds_2g[[i+1]] <- combined
      }
    }
    inds_2g <- inds_2g[!is.na(inds_2g)]
  }
  message("Finished sampling CpG clusters.")
  message("Five-number summary of number of CpGs in VMRs: ")
  message(round(quantile(lengths(inds_2g))), domain = NA, appendLF = TRUE)
  
  # Infer indices of null segments in between VMRs
  inds_1g <- list()
  inds_1g[[1]] <- 1:(inds_2g[[1]][1]-1); 
  for (i in 1:length(inds_2g)) {
    h <- inds_2g[[i]][length(inds_2g[[i]])]+1
    r <- ifelse(i < length(inds_2g), yes = inds_2g[[i+1]][1]-1, no = length(total))
    cat("h - r = ", h - r, "; h = ", h, "; r = ", r, "\n")
    if (h <= r) inds_1g[[i+1]] <- h:r 
  }
  message("Five-number summary of number of null CpGs between neighboring VMRs: ")
  message(round(quantile(lengths(inds_1g))), domain = NA, appendLF = TRUE)
  message("Finished sampling CpG clusters.")
  
  # Sample prevalences for VMRs
  pi_prob <- choose(NP, 1:(NP-1))/sum(choose(NP, 1:(NP-1)))
  pi_space <- (1:(NP-1))*1/NP
  pi1s <- sample(pi_space, size = NV, replace = TRUE, prob = pi_prob)
  
  # Pre-computation
  par_u <- vmrseq:::.priorParams(med_cov = median(total), type = "u")
  par_m <- vmrseq:::.priorParams(med_cov = median(total), type = "m")
  pars <- list("u" = par_u, "m" = par_m)
  
  # Generate hidden states
  state_seqs_1g <- lapply(inds_1g, function(ix) state[ix]) 
  state_seqs_2g <- lapply(inds_2g, function(ix) rep(1, length(ix))) # ps: hidden states in VMR are all (1,0)
  

  # Modify real-data M matrix to model-simulated data
  for (i in 1:length(inds_2g)) {
    ix <- inds_2g[[i]]
    miss_mat <- ifelse(M_mat[ix, ] >= 0, 1, NA)
    M_mat[ix, ] <- .sampScMeth2Grp(
      state_seq = state_seqs_2g[[i]],
      pi1 = pi1s[i],
      miss_mat, 0, pars
    )
    cat(i, " ")
  }
  for (i in 1:length(inds_1g)) {
    ix <- inds_1g[[i]]
    if (!is.null(ix)) {
      miss_mat <- ifelse(M_mat[ix, ] >= 0, 1, NA)
      if (is.null(dim(miss_mat))) miss_mat <- matrix(miss_mat, ncol = N)
      M_mat[ix, ] <- .sampScMeth1Grp(
        state_seq = state_seqs_1g[[i]],
        miss_mat, 0, pars
      )
    }
    cat(i, " ")
  }
  message("Finished filling in single-cell methyl values for VMRs.")

  # Useful metadata for CpG sites
  gr$is_vml <- rep(FALSE, length(gr))
  gr$vmr_name <- rep(NA, length(gr))
  gr$pi1 <- rep(NA, length(gr))
  for (i in 1:length(inds_2g)) {
    ix <- inds_2g[[i]]
    gr$is_vml[ix] <- TRUE
    gr$vmr_name[ix] <- i
    gr$pi1[ix] <- pi1s[i]
  }
  gr$meth <- rowSums(M_mat, na.rm = TRUE)
  gr$total <- rowSums(!is.na(M_mat))
  gr$mf <-   gr$meth / gr$total
  message("Finished adding metadata for CpGs.")

  # Initialize SummrizedExperiment object of the pseudo chr
  write_dir <- paste0(out_dir, "simChr_",
                      subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
                      NV, "VMRs_seed", seed)
  saveHDF5SummarizedExperiment(
    x = SummarizedExperiment(
      assays = list("M_mat" = M_mat),
      rowRanges = gr
    ),
    dir = write_dir,
    replace = TRUE
  )

  message(paste0("Finished saving simulated chromosome in a HDF5SummarizedExperiment object to:\n '",
                 write_dir, "'"))
  return(write_dir)
}


# generateNullPseudoChr <- function(
#     N, # total number of cells
#     NV, # number of VMRs
#     gamma = 0.9, # methylation info missing rate
#     sigma = 0, # noise level
#     seed = 2022 # random seed
# ) {
#   
#   set.seed(seed)
#   
#   # Pre-computation
#   par_u <- vmrseq:::.priorParams(med_cov = N*(1-gamma), type = "u")
#   par_m <- vmrseq:::.priorParams(med_cov = N*(1-gamma), type = "m")
#   pars <- list("u" = par_u, "m" = par_m)
#   
#   # Sample number of CpGs of VMRs and null regions (region between VMRs)
#   nc_1g <- rnbinom(NV + 1, size = 50, prob = 1/20)
# 
#   # Infer list of indices for each VMR and null region
#   inds_1g <- lapply(nc_1g, seq_len)
#   flag <- inds_1g[[1]][length(inds_1g[[1]])]
#   for (i in 1:NV) {
#     inds_1g[[i+1]] <- inds_1g[[i+1]] + flag
#     flag <- inds_1g[[i+1]][length(inds_1g[[i+1]])]
#   }
#   
#   # Generate genomic positions by pulling a random segment from real data
#   NC <- sum(nc_1g) # = total number of CpGs in pseudo chr
#   ind_head <- sample(x = 1:(length(dat0$pos)-NC+1), size = 1) # = pseudo-chr start index in read data
#   pos <- dat0$pos[ind_head:(ind_head+NC-1)]
#   pos <- pos - pos[1] + 1 # make pos start from 1
#   state <- dat0$state[ind_head:(ind_head+NC-1)] # reference hidden states from real data
#   
#   # Generate hidden states
#   state_seqs_1g <- lapply(inds_1g, function(ix) state[ix]) 
# 
#   # Generate MF matrix for regions
#   MFs_1g <- bplapply(1:(NV+1), function(i) .sampScMeth1Grp(state_seqs_1g[[i]], N, gamma, sigma, pars))
# 
#   # Useful metadata for CpG sites
#   is_vml <- rep(FALSE, NC)
#   vmr_name <- rep(NA, NC)
#   pi1 <- rep(NA, NC)
# 
#   # Initialize SummrizedExperiment object of the pseudo chr
#   SE <- SummarizedExperiment(
#     rowRanges = GRanges(
#       seqnames = "pseudo", 
#       ranges = IRanges(start = pos, end = pos), 
#       is_vml = is_vml, 
#       vmr_name = vmr_name, 
#       pi1 = pi1
#     ),
#     assays = list(
#       "MF" = do.call(rbind, lapply(1:(NV+1), function(i) rbind(MFs_1g[[i]])))
#     )
#   )
#   
#   return(SE)
# }




