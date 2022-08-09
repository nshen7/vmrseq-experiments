.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
# here::i_am("code/sim_studies/helper_functions/helper_functions_generate_chr.R")
source("code/sim_studies/benchmark_sim_chr/helper_functions_generate_region.R")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(SummarizedExperiment))

# register(SnowParam(workers = 6, type = "FORK"))

# Import CpG distances in chr 1 from a real dataset
dat0 <- fread("data/processed/summarized_liu2021/summarized_subtype_MSN-D1_Khdrbs3_540cells.csv.gz") %>%
  filter(chr == "chr1") %>%
  mutate(state = as.numeric(cell_meth/cell_cov >= 0.5)) %>%
  dplyr::select(c(pos, state))

# ==== Main function for generating  object of a pseudo chromosome ====

generatePseudoChr <- function(
    N, # total number of cells
    NV, # number of VMRs
    NP, # total number of subpopulations 
    seed = 2022, # random seed
    nb_r_1g = 4, # size in rnbinom function for number of contiguous null CpGs 
    nb_p_1g = 1/250, # prob in rnbinom function for number of contiguous null CpGs 
    nb_r_2g = 1, # size in rnbinom function for number of contiguous non-null CpGs 
    nb_p_2g = 1/50, # prob in rnbinom function for number of contiguous non-null CpGs 
    gamma = 0.9, # methylation info missing rate
    sigma = 0 # noise level
) {
  
  set.seed(seed)
  
  # Pre-computation
  par_u <- vmrseq:::.priorParams(med_cov = N*(1-gamma), type = "u")
  par_m <- vmrseq:::.priorParams(med_cov = N*(1-gamma), type = "m")
  pars <- list("u" = par_u, "m" = par_m)

  # Sample prevalences for VMRs
  pi_prob <- choose(NP, 1:(NP-1))/sum(choose(NP, 1:(NP-1)))
  pi_space <- (1:(NP-1))*1/NP
  pi1s <- sample(pi_space, size = NV, replace = TRUE, prob = pi_prob)

  # Sample number of CpGs of VMRs and null regions (region between VMRs)
  nc_1g <- rnbinom(NV + 1, size = nb_r_1g, prob = nb_p_1g)
  nc_2g <- rnbinom(NV, size = nb_r_2g, prob = nb_p_2g) %>% pmax(5)

  # Infer list of indices for each VMR and null region
  inds_1g <- lapply(nc_1g, seq_len)
  inds_2g <- lapply(nc_2g, seq_len)
  flag <- inds_1g[[1]][length(inds_1g[[1]])]
  for (i in 1:NV) {
    inds_2g[[i]] <- inds_2g[[i]] + flag
    flag <- inds_2g[[i]][length(inds_2g[[i]])]
    inds_1g[[i+1]] <- inds_1g[[i+1]] + flag
    flag <- inds_1g[[i+1]][length(inds_1g[[i+1]])]
  }
    
  # Generate genomic positions by pulling a random segment from real data
  NC <- sum(nc_1g) + sum(nc_2g) # = total number of CpGs in pseudo chr
  ind_head <- sample(x = 1:(length(dat0$pos)-NC+1), size = 1) # = pseudo-chr start index in read data
  pos <- dat0$pos[ind_head:(ind_head+NC-1)]
  pos <- pos - pos[1] + 1 # make pos start from 1
  state <- dat0$state[ind_head:(ind_head+NC-1)] # reference hidden states from real data
  
  # Generate hidden states
  state_seqs_1g <- lapply(inds_1g, function(ix) state[ix]) 
  state_seqs_2g <- lapply(inds_2g, function(ix) rep(1, length(ix))) # ps: hidden states in VMR are all (1,0)
  
  # Generate MF matrix for regions
  MFs_1g <- bplapply(1:(NV+1), function(i) .sampScMeth1Grp(state_seqs_1g[[i]], N, gamma, sigma, pars))
  MFs_2g <- bplapply(1:NV, function(i) .sampScMeth2Grp(state_seqs_2g[[i]], pi1s[i], N, gamma, sigma, pars))
  
  # Useful metadata for CpG sites
  is_vml <- rep(FALSE, NC)
  vmr_name <- rep(NA, NC)
  pi1 <- rep(NA, NC)
  for (i in 1:NV) {
    is_vml[inds_2g[[i]]] <- TRUE
    vmr_name[inds_2g[[i]]] <- i
    pi1[inds_2g[[i]]] <- pi1s[i]
  }
  
  # Initialize SummrizedExperiment object of the pseudo chr
  SE <- SummarizedExperiment(
    rowRanges = GRanges(
      seqnames = "pseudo", 
      ranges = IRanges(start = pos, end = pos), 
      is_vml = is_vml, 
      vmr_name = vmr_name, 
      pi1 = pi1
    ),
    assays = list(
      "MF" = rbind(MFs_1g[[1]], do.call(rbind, lapply(1:NV, function(i) rbind(MFs_2g[[i]], MFs_1g[[i+1]]))))
    )
  )
  
  return(SE)
}


generateNullPseudoChr <- function(
    N, # total number of cells
    NV, # number of VMRs
    gamma = 0.9, # methylation info missing rate
    sigma = 0, # noise level
    seed = 2022 # random seed
) {
  
  set.seed(seed)
  
  # Pre-computation
  par_u <- vmrseq:::.priorParams(med_cov = N*(1-gamma), type = "u")
  par_m <- vmrseq:::.priorParams(med_cov = N*(1-gamma), type = "m")
  pars <- list("u" = par_u, "m" = par_m)
  
  # Sample number of CpGs of VMRs and null regions (region between VMRs)
  nc_1g <- rnbinom(NV + 1, size = 50, prob = 1/20)

  # Infer list of indices for each VMR and null region
  inds_1g <- lapply(nc_1g, seq_len)
  flag <- inds_1g[[1]][length(inds_1g[[1]])]
  for (i in 1:NV) {
    inds_1g[[i+1]] <- inds_1g[[i+1]] + flag
    flag <- inds_1g[[i+1]][length(inds_1g[[i+1]])]
  }
  
  # Generate genomic positions by pulling a random segment from real data
  NC <- sum(nc_1g) # = total number of CpGs in pseudo chr
  ind_head <- sample(x = 1:(length(dat0$pos)-NC+1), size = 1) # = pseudo-chr start index in read data
  pos <- dat0$pos[ind_head:(ind_head+NC-1)]
  pos <- pos - pos[1] + 1 # make pos start from 1
  state <- dat0$state[ind_head:(ind_head+NC-1)] # reference hidden states from real data
  
  # Generate hidden states
  state_seqs_1g <- lapply(inds_1g, function(ix) state[ix]) 

  # Generate MF matrix for regions
  MFs_1g <- bplapply(1:(NV+1), function(i) .sampScMeth1Grp(state_seqs_1g[[i]], N, gamma, sigma, pars))

  # Useful metadata for CpG sites
  is_vml <- rep(FALSE, NC)
  vmr_name <- rep(NA, NC)
  pi1 <- rep(NA, NC)

  # Initialize SummrizedExperiment object of the pseudo chr
  SE <- SummarizedExperiment(
    rowRanges = GRanges(
      seqnames = "pseudo", 
      ranges = IRanges(start = pos, end = pos), 
      is_vml = is_vml, 
      vmr_name = vmr_name, 
      pi1 = pi1
    ),
    assays = list(
      "MF" = do.call(rbind, lapply(1:(NV+1), function(i) rbind(MFs_1g[[i]])))
    )
  )
  
  return(SE)
}




