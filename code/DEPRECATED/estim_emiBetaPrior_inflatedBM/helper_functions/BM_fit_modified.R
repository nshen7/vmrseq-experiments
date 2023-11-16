## install if not
# devtools::install_github("https://github.com/palmerimatthew/BetaMixture")
suppressPackageStartupMessages(library(BetaMixture))

# ==== utils ====
#Expectation step ----
E_Step <- function(x, Alphas, Betas) {
  K = length(Alphas)
  responsibility = matrix(nrow = K, ncol = length(x))
  #matrix of responsibility numerators (KxN)
  for (i in 1:K) {
    responsibility[i,] = dbeta(x, Alphas[i], Betas[i])
  }
  respons_denom = colSums(responsibility)
  #responsibility matrix (NxK)
  responsibility = t(responsibility)/respons_denom
  Mixture_Params = colSums(responsibility)/sum(responsibility)
  return(Mixture_Params)
}

#Maximization step ----
M_Step <- function(x, Mixture_Params, Old_Component_Params) {
  stopifnot("Wrong number of parameters in `Component_Prams." = length(Old_Component_Params) == 4)
  #likelihood equation for numerical MLEs
  likelihood <- function(theta, x) {
    if (any(theta <= 0)) {return(NA)}
    alphas = theta[seq(1, length(theta), by=2)]
    betas = theta[seq(2, length(theta), by=2)]
    return(BM_Density(x, Mixture_Params, alphas, betas, Log = T))
  }
  # ## ineqA and ineqB for inequality constraints A %*% theta + B > 0.
  # A <- matrix(c(0,1,0,0), 1, 4)
  # B <- -1.0001
  New_Component_Params = maxLik::maxLik(likelihood, 
                                        start = Old_Component_Params, 
                                        # constraints = list(ineqA = A, ineqB = B),
                                        x = x)
  return(New_Component_Params$estimate)
}

# # Method of Moments initial values for component parameters ----
# MoM_Calculation <- function(sorted_x, K, len) {
#   cut_points = seq(1, len, length.out = K+1)
#   MoM_Estimates = c()
#   for(i in 1:K) {
#     temp_data = sorted_x[ceiling(cut_points[i]):floor(cut_points[i+1])]
#     temp <- MoM_Beta(temp_data)
#     alpha = temp$alpha
#     beta = temp$beta
#     MoM_Estimates = c(MoM_Estimates, c(alpha, beta))
#   }
#   return(MoM_Estimates)
# }

MoM_Beta <- function(data) {
  mean = mean(data)
  var = var(data)
  intermediate = mean*(1-mean)/var - 1
  alpha = mean*intermediate
  beta = (1-mean)*intermediate
  return(list(alpha = alpha, beta = beta))
}

#Testing ----
Test_BM_Fit <- function() {
  #Precondition tests for BM_Fit ----
  print("BM_Fit tests not implemented yet")
  
  #MoM_Calculation tests ----
  print("MoM_Calculation tests not implemented yet")
}


###############################################################
#### Modified the initial parameter values to MoM estimate ####
#### from methylated and unmethylated subpopulation!!!!!!! ####
###############################################################

BM_Fit <- function(x, threshold = 0.0001, seed = NA) {

  stopifnot("`x` should be > 0 and < 1." = all(x>0 & x<1))
  #setup ----
  if (!is.na(seed)) set.seed(seed)
  
  K <- 2
  
  x_u <- x[x < 0.5]; x_m <- x[x > 0.5]
  w_u <- length(x_u) / length(x); w_m <- 1-w_u
  
  Mix_Params = c(w_u, w_m)
  Component_Params = unname(c(unlist(MoM_Beta(x_u)), unlist(MoM_Beta(x_m))))
  Component_Params = M_Step(x, Mix_Params, Component_Params)
  Alphas = Component_Params[seq(1, length(Component_Params), by=2)]
  Betas = Component_Params[seq(2, length(Component_Params), by=2)]
  current_log_likelihood = sum(BM_Density(x, Mix_Params, Alphas, Betas, T))
  previous_log_likelihood = current_log_likelihood - 100
  #record of iterations
  iteration_results = matrix(nrow = 100, ncol = 2+3*K)
  iteration = 1
  iteration_results[1,] = c(iteration, Mix_Params, Alphas, Betas, current_log_likelihood)
  #iterative process ----
  while (abs(current_log_likelihood - previous_log_likelihood) > threshold) {
    if(iteration >= 100) break
    iteration = iteration + 1
    #E step
    Mix_Params = E_Step(x, Alphas, Betas)
    #M step
    Component_Params = M_Step(x, Mix_Params, Component_Params)
    Alphas = Component_Params[seq(1, length(Component_Params), by=2)]
    Betas = Component_Params[seq(2, length(Component_Params), by=2)]
    previous_log_likelihood = current_log_likelihood
    current_log_likelihood = sum(BM_Density(x, Mix_Params, Alphas, Betas, T))
    iteration_results[iteration,] = c(iteration, Mix_Params, Alphas, Betas, current_log_likelihood)
    cat("current_log_likelihood: ", current_log_likelihood, "previous_log_likelihood: ", previous_log_likelihood, "\n")
  }
  #preparing stuff for return statements ----
  iteration_results <- as.data.frame(iteration_results)
  iteration_results <- dplyr::filter(iteration_results, !is.na(V1))
  #making column names
  columns <- "Iteration"
  for (i in 1:K) {columns = c(columns, paste0('Mixture_Param_', i))}
  for (i in 1:K) {columns = c(columns, paste0('Alpha_', i))}
  for (i in 1:K) {columns = c(columns, paste0('Beta_', i))}
  columns <- c(columns, "Log_Likelihood")
  colnames(iteration_results) <- columns
  #returning results
  return(list(Mix_Params = Mix_Params, Alpha = Alphas, Beta = Betas, Iterations = iteration_results, LogLik = current_log_likelihood))
}
