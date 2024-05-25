source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")

# function for computing nearst neighbor score (reference: scbs paper)
nnScore <- function(d_mat, true_clust, k, theta) {
  
  if(ncol(d_mat) != nrow(d_mat)) stop("`d_mat` is not square")
  
  nnRow <- function(i) {
    d_vec <- d_mat[i,]
    idx_nn <- order(d_vec)[2:(k+1)] # first element is self so removed
    
    clust_nn <- true_clust[idx_nn]
    clust_self <- true_clust[i]
    
    return(sum(clust_nn == clust_self) / k) # return the proportion of correctly assigned NNs
  }
  
  N <- nrow(d_mat)
  props <- do.call(c, lapply(1:N, nnRow))
  return(sum(props > theta) / N)
}
