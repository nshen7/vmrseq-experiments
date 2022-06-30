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
here::i_am("code/package_functions/helper_functions.R")
source(here::here("code/package_functions/helper_functions.R"))

#### Example ####

## pre-computing
K <- 15; lambda <- 20; N <- 100; gamma <- 0.93; sigma <- 0.05; seed <- 1000
# inits <- c(0.25, 0.5, 0.75)
inits <- c(0.1, 0.5, 0.9)
REFARRAY <- .calRefArray(max_size = N)
CHOICEARRAY <- .calChoiceArray(REFARRAY)
list <- .calMethArray(par_u = .priorParams(med_cov = N*(1-gamma), type = "u"),
                      par_m = .priorParams(med_cov = N*(1-gamma), type = "m"),
                      max_size = N)
METHARRAY <- list$METHARRAY; UNMETHARRAY <- list$UNMETHARRAY
pos <- .sampGenomCoord(K, lambda, seed)

## 1 grouping
seed <- 1
(state_seq_1g <- .sampHiddState1Grp(pos, seed))
MF_1g <- .sampScMeth1Grp(state_seq_1g, N, gamma, sigma, seed)
(meth_reads <- rowSums(MF_1g, na.rm = T))
(total_reads <- rowSums(!is.na(MF_1g)))
round(meth_reads / total_reads, 2)
# plot(pos, meth_reads/total_reads, ylim = c(0,1)); sum(meth_reads/total_reads < 1)
# plot(pos, values(cells.se)$cell_MF[201:300], ylim = c(0,1)); sum(values(cells.se)$cell_MF[201:300] < 1)

.Viterbi1Grp(pos, total_reads, meth_reads, tp = NULL, METHARRAY, UNMETHARRAY)
system.time(.prevOptimMultiInit(pos, total_reads, meth_reads, inits,  backtrack = T,
                    CHOICEARRAY = CHOICEARRAY, 
                    METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY))

## 2 grouping
pi1 <- 0.3; seed <- 1
(state_seq_2g <- .sampHiddState2Grp(pos, seed))
MF_2g <- .sampScMeth2Grp(state_seq_2g, pi1, N, gamma, sigma, seed)
(meth_reads <- rowSums(MF_2g, na.rm = T))
(total_reads <- rowSums(!is.na(MF_2g)))
round(meth_reads / total_reads, 2)
plot(pos, meth_reads/total_reads, ylim = c(0,1))

res_1g <- .Viterbi1Grp(pos, total_reads, meth_reads, tp = NULL, METHARRAY, UNMETHARRAY)
res_1g[K, 2]
system.time(res_2g <- .prevOptimMultiInit(pos, total_reads, meth_reads, inits, backtrack = T,
                              CHOICEARRAY = CHOICEARRAY, 
                              METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY))
res_2g$loglik; res_2g$optim_pi_1; rowSums(res_2g$vit_path)
## seems that backtracking works well when K is small, not so well when K is large


# ==== utils ====

.normTo1 <- function(x) x / sum(x)

### Sample `p` of length n from beta-mixture priors, where the prior distributions vary depending on n.
.sampBeta <- function(n, med_cov, state) {
  if (!state) {
    par_u <- .priorParams(med_cov = med_cov, type = "u")
    return(rBEZI(n, mu = par_u['mu'], sigma = par_u['sigma'], nu = par_u['nu']))
    # return(replicate(n, ifelse(runif(1)<=par_u[1],
    #                            yes = 0, 
    #                            no = rbeta(1, par_u[2], par_u[3])))
    # )
  } else {
    par_m <- .priorParams(med_cov = med_cov, type = "m")
    return(rBE(n, mu = par_m['mu'], sigma = par_m['sigma']))
    # return(replicate(n, ifelse(runif(1)<=par_m[1],
    #                            yes = 1,
    #                            no = rbeta(1, par_m[2], par_m[3])))
    # )
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
  state_seq[1] <- sample(0:1, 1, prob = c(0.2, 0.8))
  
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
  nois_seq <- rnorm(n = K, 0, sigma)
  p0_seq <- map_dbl(state_seq, ~.sampBeta(n = 1, med_cov = N*(1-gamma), state = .x))
  p_mat <- matrix(replicate(N, (p0_seq + nois_seq) %>% pmax(0) %>% pmin(1)), ncol = N)
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
  state_seq[1] <- sample(0:2, 1, prob = c(0.2^2, 0.2*0.8, 0.8^2) %>% .normTo1())
  
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









