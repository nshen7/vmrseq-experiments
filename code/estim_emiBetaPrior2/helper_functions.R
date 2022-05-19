library(data.table)

# ==== complete function (estimate params from MF vector) ====

### `population` == "m" means estimate prior parameters for methylated population, or "u" for unmethylated population
calBestParam <- function(mf_vec, 
                         population,
                         n_bins = 200, 
                         param1_space = 1:20, 
                         param2_space = seq(100,2000,10), 
                         w1_space = seq(0.6,1,0.01),
                         mc.cores){
  
  dens <- calDensity(mf_vec, n_bins)

  if (population == "u"){
    dens[x >= 0.5, "y"] <- 0
    dens$y <- dens$y / sum(dens$y, na.rm = T)
  } else if (population == "m") {
    dens[x <= 0.5, "y"] <- 0
    dens$y <- dens$y / sum(dens$y, na.rm = T)
  } else {
    stop("population should be either \"m\" or \"u\".")
  }
  
  return(calBestChoice(dens, population, 
                       param1_space, param2_space,  w1_space,
                       mc.cores))
}



# ==== partial function (estimate params from density) ====
calBestParamFromDensity <- function(dens, # density estimated from `calDensity`
                                    population,
                                    param1_space = 1:20, 
                                    param2_space = seq(100,2000,10), 
                                    w1_space = seq(0.6,1,0.01),
                                    mc.cores){
  n_bins <- nrow(dens) - 1
  dens <- as.data.table(dens)
  colnames(dens) <- c("x", "y")

  if (population == "u"){
    dens[x >= 0.5, "y"] <- 0
    dens$y <- dens$y / sum(dens$y, na.rm = T)
  } else if (population == "m") {
    dens[x <= 0.5, "y"] <- 0
    dens$y <- dens$y / sum(dens$y, na.rm = T)
  } else {
    stop("population should be either \"m\" or \"u\".")
  }
  
  return(calBestChoice(dens, population,
                       param1_space, param2_space,  w1_space,
                       mc.cores))
}

# === TODO: plotting functions ==== 


# ==== utilities ====


calDensity <- function(vec, n_bins = 200) { 
  
  df <- data.table("x" = seq(0, 1, 1/n_bins), "y" = 0)
  tab <- table(round(vec*n_bins)/n_bins) # methylation fraction rounded to the nearest 0.5% if n_bins=200
  
  df[, "y"] <- sapply(1:nrow(df), function(i) tab[as.character(df[i, "x"])]/sum(tab))
  # df[is.na(y), "y"] <- 0 ### changed
  return(df)
}


calL1Diff <- function(dens, population, w1, param1, param2, n_bins){
  
  if (population == "u"){
    beta <- w1*dbeta(seq(0, 1, 1/n_bins), 1, param1) + (1-w1) * dbeta(seq(0, 1, 1/n_bins), 1, param2)
  } else if (population == "m") {
    beta <- w1*dbeta(seq(0, 1, 1/n_bins), param1, 1) + (1-w1) * dbeta(seq(0, 1, 1/n_bins), param2, 1)
  }
  
  beta <- beta/sum(beta)
  return(sum(abs(beta-dens$y), na.rm = T))
}


calBestChoice <- function(dens, population, 
                          param1_space, param2_space,  w1_space,
                          mc.cores) {
  n_bins <- nrow(dens) - 1
  param_space <- expand.grid(w1_space, param1_space, param2_space)
  colnames(param_space) <- c("w1", "param1", "param2")
  param_space$l1_diff <- mclapply(1:nrow(param_space), 
                                  function(i) calL1Diff(dens, population,
                                                        param_space$w1[i], 
                                                        param_space$param1[i], 
                                                        param_space$param2[i], 
                                                        n_bins),
                                  mc.cores = mc.cores
  )
  best_choice <- param_space[which.min(param_space$l1_diff),]
  
  if (population == "u"){
    cat(paste0(
      "Fitted mixture distribution density for unmethylated population: ", best_choice$w1 ,
      " * dbeta(x, 1, ", best_choice$param1, ") + ", round(1-best_choice$w1, 2), 
      " * dbeta(x, 1, ", best_choice$param2, ") \n",
      "L1 difference between discrete empirical distribution and fitted mixture:", best_choice$l1_diff, "\n")
    )
  } else if (population == "m") {
    cat(paste0(
      "Fitted mixture distribution density for methylated population: ", best_choice$w1 ,
      " * dbeta(x, ", best_choice$param1, ", 1) + ", round(1-best_choice$w1, 2), 
      " * dbeta(x, ", best_choice$param2, ", 1) \n",
      "L1 difference between discrete empirical distribution and fitted mixture:", best_choice$l1_diff, "\n")
    )
  }
  
  
  # return(list(best_param1 = best_choice$param1, 
  #             best_param2 = best_choice$param2, 
  #             best_w = best_choice$param2, 
  #             l1_diff = best_choice$l1_diff))
  return(best_choice)
}
