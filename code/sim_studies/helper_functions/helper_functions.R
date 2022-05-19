# - Variables of interest:
#   - Number of CpGs in region $K$
#   - Number of cells $N$
#   - Level of CpG density
#     - Rich: CpG-CpG distance ~ $\text{Poisson}(\lambda = 10)$
#     - Moderate:  CpG-CpG distance ~ $\text{Poisson}(\lambda = 250)$
#     - Poor: CpG-CpG distance ~ $\text{Poisson}(\lambda = 1,000)$
#   - Methylation info missing rate $\gamma$ (hence average across-cell coverage should be $N\gamma$ )
#   - Noise level $\sigma$
#   - Prevalence of the methylated grouping $\pi_1$ (only for 2-grouping case)

suppressPackageStartupMessages(library(tidyverse))
source(here::here("code/package_functions/helper_functions.R"))

#### Example ####

## pre-computing
K <- 10; lambda <- 10; N <- 300; gamma <- 0.8; sigma <- 0.2; seed <- 2001
pi1 <- 0.3; inits <- c(0.25, 0.5, 0.75)
pos <- .sampGenomCoord(K, lambda, seed)
REFARRAY <- .calRefArray(max_size = N)
CHOICEARRAY <- .calChoiceArray(REFARRAY)
### For 1-grouping
list <- .calMethArray(par_u = .priorParams(N, type = "u"),
                      par_m = .priorParams(N, type = "m"),
                      REFARRAY)
METHARRAY <- list$METHARRAY; UNMETHARRAY <- list$UNMETHARRAY
### For 2-grouping
list <- map(inits, ~ .calMethArray(par_u = .priorParams(round(N*(1-.x)), type = "u"),
                                   par_m = .priorParams(round(N*.x), type = "m"),
                                   REFARRAY))
METHLIST <- map(list, ~.x$METHARRAY); UNMETHLIST <- map(list, ~.x$UNMETHARRAY) # to store the beta-binomial prob matrix (no log)

## 1 grouping
(state_seq_1g <- .sampHiddState1Grp(pos, seed))
MF_1g <- .sampScMeth1Grp(state_seq_1g, N, gamma, sigma, seed)
(meth_reads <- rowSums(MF_1g, na.rm = T))
(total_reads <- rowSums(!is.na(MF_1g)))
round(meth_reads / total_reads, 2)
.Viterbi1Grp(pos, total_reads, meth_reads, tp = NULL, METHARRAY, UNMETHARRAY)
.prevOptimMultiInit(pos, total_reads, meth_reads, inits, 
                    CHOICEARRAY = CHOICEARRAY, METHLIST = METHLIST, UNMETHLIST = UNMETHLIST)

## 2 grouping
(state_seq_2g <- .sampHiddState2Grp(pos, seed))
MF_2g <- .sampScMeth2Grp(state_seq_2g, pi1, N, gamma, sigma, seed)
(meth_reads <- rowSums(MF_2g, na.rm = T))
(total_reads <- rowSums(!is.na(MF_2g)))
round(meth_reads / total_reads, 2)
.Viterbi1Grp(pos, total_reads, meth_reads, tp = NULL, METHARRAY, UNMETHARRAY)
.prevOptimMultiInit(pos, total_reads, meth_reads, inits, 
                    CHOICEARRAY = CHOICEARRAY, METHLIST = METHLIST, UNMETHLIST = UNMETHLIST)


# ==== utils ====

.normTo1 <- function(x) x / sum(x)

### Sample p of length n from beta-mixture priors, where the prior distributions vary depending on n.
.sampBeta <- function(n, state) {
  if (!state) {
    par_u <- .priorParams(N = n, type = "u")
    return(replicate(n, ifelse(runif(1)<=par_u[1], 
                               yes = rbeta(1, 1, par_u[2]), 
                               no = rbeta(1, 1, par_u[3])))
    )
  } else {
    par_m <- .priorParams(N = n, type = "m")
    return(replicate(n, ifelse(runif(1)<=par_m[1], 
                               yes = rbeta(1, par_m[2], 1), 
                               no = rbeta(1, par_m[3], 1)))
    )
  }
}

### Sample missing status (of N cells) for each CpG site
.sampMiss1Cell <- function(N, gamma) sample(c(NA, 1), N, prob = c(gamma, 1-gamma), replace = T)


# ==== Functions for sampling NULL regions ====

### Sample genomic coordinates
.sampGenomCoord <- function(K, lambda, seed) {
  set.seed(seed)
  return(cumsum(rpois(K, lambda = lambda)))
}

### Sample undeerlying methylation state for 1-grouping case
.sampHiddState1Grp <- function(pos, seed, tp_mat = .loadTransitProbs(pos)) {
  K <- length(pos)
  stopifnot("Dimension of `tp_mat` not matched to number of CpGs `K`." = nrow(tp_mat) == K-1)
  
  ## Initialization
  set.seed(seed)
  state_seq <- rep(-1, K) 
  
  ## Sample initial state 
  state_seq[1] <- sample(0:1, 1, prob = c(0.23, 0.77))
  
  ## Sample subsequent states based on transition probability
  for (i in 2:K) {
    probs <- tp_mat[i-1, state_seq[i-1] + c(1,3)] %>% unlist()
    state_seq[i] <- sample(0:1, 1, prob = probs)
  }
  return(state_seq)
}

### Sample methylation level for individual cells
.sampScMeth1Grp <- function(state_seq, N, gamma, sigma, seed) {
  K <- length(state_seq)
  stopifnot("Number of CpGs `K` not equal to length of state sequence." = K == length(state_seq))
  stopifnot("Gamma should be between 0 and 1." = gamma > 0 & gamma < 1)
  set.seed(seed)
  
  ## sample missing status
  miss_mat <- t(replicate(K, .sampMiss1Cell(N, gamma)))
  ## sample probability of being methylated
  nois_mat <- matrix(data = runif(K*N, -sigma, sigma), nrow = K, ncol = N)
  p0_mat <- do.call(rbind, map(state_seq, ~.sampBeta(n = N, state = .x))) %>% as.matrix()
  p_mat <- (p0_mat + nois_mat) %>% pmax(0) %>% pmin(1)
  ## sample methylation levels
  MF <- apply(p_mat, c(1, 2), function(p) rbernoulli(1, p = p) %>% as.integer) * miss_mat
  
  return(MF)
}



# ==== Functions for sampling NON-NULL regions ====

### Sample underlying methylation state for 2-grouping case
.sampHiddState2Grp <- function(pos, seed, tp_mat = .loadTransitProbs(pos)) { 
  K <- length(pos)
  stopifnot("Dimension of `tp_mat` not matched to number of CpGs `K`." = nrow(tp_mat) == K-1)
  
  ## Initialization
  set.seed(seed)
  state_seq <- rep(-1, K) ## states in integer representation: `0`->(0,0), `1`->(1,0), `2`->(1,1)
  
  ## Sample initial state 
  state_seq[1] <- sample(0:2, 1, prob = c(0.23^2, 0.23*0.77, 0.77^2) %>% .normTo1())
  
  ## Sample subsequent states based on transition probability
  for (i in 2:K) {
    ## compute transition prob dist to next state
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
    ## sample next state
    state_seq[i] <- sample(0:2, 1, prob = probs)
  }
  return(state_seq)
}


.sampScMeth2Grp <- function(state_seq, pi1, N, gamma, sigma, seed, tp_mat = .loadTransitProbs(pos)) {
  state_vec <- do.call(rbind, map(state_seq, .translateState2Grp))
  MF_u <- .sampScMeth1Grp(state_vec[,1], ceiling(pi1*N), gamma, sigma, seed)
  MF_m <- .sampScMeth1Grp(state_vec[,2], N-ceiling(pi1*N), gamma, sigma, seed)
  return(cbind(MF_u, MF_m))
}









