source("code/SETPATHS.R")
library(SummarizedExperiment)
library(HDF5Array)
source("code/case_studies/util_functions/6nnScore.R")

read_dir <- "data/interim/case_studies/luo2017mice_full_sens_analysis/result_summary/minN"
write_dir <- "data/interim/case_studies/luo2017mice_full_sens_analysis/result_summary/minN"
if (!file.exists(write_dir)) dir.create(write_dir)
md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")

# ---- utils ----

computeDissimilarityMatrix <- function(se, dissim_metric = "manhattan") {
  MF <- as.matrix(assays(se)$M/assays(se)$Cov) 
  d_mat <- cluster::daisy(t(MF), metric = dissim_metric, stand = FALSE) %>% as.matrix()
  return(d_mat)
}

computeNNScoreFromSE <- function(minN, k, theta) {

  vmrs.se <- loadHDF5SummarizedExperiment(here(read_dir, paste0('minN', minN), "vmrseq_regionSummary_vmrs"))
  d_mat_vmrs <- computeDissimilarityMatrix(vmrs.se)
  fwrite(d_mat_vmrs, here(write_dir, paste0("dissimilarity_matrix_regional_methyl_minN", minN, "_vseq.txt.gz")), 
         col.names = F, row.names = F, quote = F)
  
  crs.se <- loadHDF5SummarizedExperiment(here(read_dir, paste0('minN', minN), "vmrseq_regionSummary_crs"))
  d_mat_crs <- computeDissimilarityMatrix(crs.se)
  fwrite(d_mat_crs, here(write_dir, paste0("dissimilarity_matrix_regional_methyl_minN", minN, "_vseq_cr.txt.gz")), 
         col.names = F, row.names = F, quote = F)
  
  NNScore_vmrs_broad <- nnScore(d_mat_vmrs, true_clust = md$Neuron_type1, k, theta)
  NNScore_crs_broad  <- nnScore(d_mat_crs, true_clust = md$Neuron_type1, k, theta)
  NNScore_vmrs_sub   <- nnScore(d_mat_vmrs, true_clust = md$Neuron.type, k, theta)
  NNScore_crs_sub    <- nnScore(d_mat_crs, true_clust = md$Neuron.type, k, theta)
  
  return(c(NNScore_vmrs_broad, NNScore_crs_broad, NNScore_vmrs_sub, NNScore_crs_sub))
}


# ---- Main ----

score.df <- data.frame(
  minN = c(3,5,10,20),
  # minN = c(0.025), 
  NNScore_vmrs_broad = 0,
  NNScore_crs_broad = 0,
  NNScore_vmrs_sub = 0,
  NNScore_crs_sub = 0
)


k <- 100
theta <- 0.7

for (i in 1:nrow(score.df)){
  score.df[i, 2:5] <- computeNNScoreFromSE(score.df$minN[i], k = k, theta = theta)
  print(i)
}

fwrite(score.df, here(write_dir, paste0("nearest_neighbor_score_k", k, "_theta", theta, ".csv")))


