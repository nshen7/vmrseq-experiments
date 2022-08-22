# - Variables of interest:
#   - Number of CpGs in region $K$
#   - Number of cells $N$
#   - Methylation info missing rate $\gamma$ (hence average across-cell coverage should be $N\gamma$ )
#   - Noise level $\sigma$
#   - Prevalence of the methylated grouping $\pi_1$ (only for 2-grouping case)
.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
# here::i_am("code/sim_studies/helper_functions/helper_functions_generate_region.R")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gamlss.dist))


# ==== utils ====

.normTo1 <- function(x) x / sum(x)

# Sample `p` of length n from beta-mixture priors, where the prior distributions vary depending on n.
# .sampBeta <- function(n, state, pars) {
#   if (!state%in%c(0,1)) stop("'state' should be either 0 or 1")
#   if (!state) {
#     par_u <- pars[['u']]
#     return(rBEZI(n, mu = par_u['mu'], sigma = par_u['sigma'], nu = par_u['nu']))
#   } else {
#     par_m <- pars[['m']]
#     return(rBE(n, mu = par_m['mu'], sigma = par_m['sigma']))
#   }
# }

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


# Sample missing status (of N cells) for 1 CpG site
.sampMiss1Cell <- function(N, gamma) {
  status <- sample(c(NA, 1), N, prob = c(gamma, 1-gamma), replace = T)
  if (all(is.na(status))) status[sample(1:N, 1)] <- 1 # At least 1 cell is covered
  return(status)
}
  
# Sample missing status (of N cells) for K CpG sites
.sampMissMat <- function(N, K, gamma) 
  do.call(rbind, replicate(K, .sampMiss1Cell(N, gamma), simplify = F))

# Sample methylation level for individual cells
.sampMethMat <- function(state_seq, miss_mat, N, sigma, pars) {
  
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

# ==== Functions for sampling NULL regions ====

# Sample undeerlying methylation state for 1-grouping case
.sampHiddState1Grp <- function(pos, tp_mat = vmrseq:::.loadTransitProbs(pos)) {
  
  K <- length(pos)
  stopifnot("Dimension of `tp_mat` not matched to number of CpGs `K`." = nrow(tp_mat) == K-1)
  
  # Initialization
  state_seq <- rep(-1, K) 
  
  # Sample initial state 
  state_seq[1] <- sample(0:1, 1, prob = c(0.2, 0.8))
  
  # Sample subsequent states based on transition probability
  for (i in 2:K) {
    probs <- tp_mat[i-1, state_seq[i-1] + c(1,3)] %>% unlist()
    state_seq[i] <- sample(0:1, 1, prob = probs)
  }
  return(state_seq)
}

# Sample methylation level for individual cells
.sampScMeth1Grp <- function(state_seq, miss_mat, sigma, pars) {
  
  N <- ncol(miss_mat)
  K <- nrow(miss_mat)
  stopifnot("Number of CpGs `K` not equal to length of state sequence." = 
              K == length(state_seq))

  mf_mat <- .sampMethMat(state_seq, miss_mat, N, sigma, pars)
  return(mf_mat * miss_mat)
}



# ==== Functions for sampling NON-NULL regions ====

# Sample underlying methylation state for 2-grouping case
.sampHiddState2Grp <- function(pos, tp_mat = .loadTransitProbs(pos)) { 
  
  K <- length(pos)
  stopifnot("Dimension of `tp_mat` not matched to number of CpGs `K`." = nrow(tp_mat) == K-1)
  
  # Initialization
  state_seq <- rep(-1, K) # states in integer representation: `0`->(0,0), `1`->(1,0), `2`->(1,1)
  
  # Sample initial state 
  state_seq[1] <- sample(0:2, 1, prob = c(0.2^2, 0.2*0.8, 0.8^2) %>% .normTo1())
  
  # Sample subsequent states based on transition probability
  for (i in 2:K) {
    # compute transition prob dist to next state
    if (state_seq[i-1] == 0) {
      probs <- c(tp_mat[i-1,1]^2, # P(0|0)P(0|0)
                 tp_mat[i-1,3]*tp_mat[i-1,1], # P(1|0)P(0|0)
                 tp_mat[i-1,3]^2) %>% .normTo1() # P(1|0)P(1|0)
    } else if (state_seq[i-1] == 1) {
      probs <- c(tp_mat[i-1,2]*tp_mat[i-1,1],
                 tp_mat[i-1,4]*tp_mat[i-1,1],
                 tp_mat[i-1,4]*tp_mat[i-1,3]) %>% .normTo1()
    } else {
      probs <- c(tp_mat[i-1,3]^2,
                 tp_mat[i-1,4]*tp_mat[i-1,3],
                 tp_mat[i-1,4]^2) %>% .normTo1()
    }
    # sample next state
    state_seq[i] <- sample(0:2, 1, prob = probs)
  }
  return(state_seq)
}


# .sampScMeth2Grp <- function(state_seq, pi1, N, gamma, sigma, pars) {
#   stopifnot("Gamma should be between 0 and 1." = gamma > 0 & gamma < 1)
#   
#   state_vec <- do.call(rbind, map(state_seq, vmrseq:::.translateState2Grp))
#   miss_mat <- .sampMissMat(N = N, K = length(state_seq), gamma)
#   ind_g1 <- 1:ceiling(pi1*N)
#   ind_g2 <- (ceiling(pi1*N)+1):N
#   mf_mat <- cbind(
#     .sampMethMat(state_seq = state_vec[,1], miss_mat = miss_mat[,ind_g1], N = length(ind_g1), sigma, pars),
#     .sampMethMat(state_seq = state_vec[,2], miss_mat = miss_mat[,ind_g2], N = length(ind_g2), sigma, pars)
#   )
#   return(mf_mat * miss_mat)
# }

.sampScMeth2Grp <- function(state_seq, pi1, miss_mat, sigma, pars) {
  
  N <- ncol(miss_mat)
  K <- nrow(miss_mat)
  
  state_vec <- do.call(rbind, map(state_seq, vmrseq:::.translateState2Grp))
  
  # ind_g1 <- 1:ceiling(pi1*N)
  # ind_g2 <- (ceiling(pi1*N)+1):N
  # mf_mat <- cbind(
  #   .sampMethMat(state_seq = state_vec[,1], miss_mat = miss_mat[,ind_g1], N = length(ind_g1), sigma, pars),
  #   .sampMethMat(state_seq = state_vec[,2], miss_mat = miss_mat[,ind_g2], N = length(ind_g2), sigma, pars)
  # )
  ind_g1 <- sample(1:N, ceiling(pi1*N)) %>% sort
  ind_g2 <- which(! 1:N %in% ind_g1)
  mf_mat <- matrix(nrow = K, ncol = N)
  mf_mat[,ind_g1] <- .sampMethMat(state_seq = state_vec[,1], miss_mat = miss_mat[,ind_g1], N = length(ind_g1), sigma, pars)
  mf_mat[,ind_g2] <- .sampMethMat(state_seq = state_vec[,2], miss_mat = miss_mat[,ind_g2], N = length(ind_g2), sigma, pars)
  return(mf_mat * miss_mat)
}












