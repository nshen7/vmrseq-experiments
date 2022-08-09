.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BiocParallel))
# suppressPackageStartupMessages(library(vmrseq))
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
devtools::load_all("../vmrseq-package/vmrseq/")
set.seed(2022)

########################################################################
######## Train transition probablities used as default in vmrseq #######
########################################################################

# Using default parameters for estimating transition probs:
max_dist_bp = 2000; buffer_bp = 5000
lags = 1:100
degree = 2; span = 0.02
print(lags)

##### Sample cells in Liu2021 data and form into list of data.frames ====

wrapper1 <- function(subtype_dir){
  
  st <- fread(subtype_dir) %>%
    filter(cell_cov >= 5) %>% # QC: across-cell coverage >= 5
    mutate(cell_MF = cell_meth / cell_cov) %>%
    mutate(state = round(cell_MF)) %>%
    group_by(chr) %>%
    mutate(lag_state = lag(state, 1))
  
  # decide the states for sites with MF = 0.5 based on its previous site
  st$state[which(st$cell_MF==0.5)] <- st$lag_state[which(st$cell_MF==0.5)]
  
  smr <- vmrseq:::.computeProb1Unit(
    df = st %>% select(chr, pos, state),
    max_dist_bp = max_dist_bp,
    buffer_bp = buffer_bp,
    lags = lags
  )

  return(smr)
}

# Information of all subtypes used for training
dir <- ("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/data/processed/summarized_liu2021/")
subtype_dirs <- paste0(dir, list.files(dir))

# Compute empirical probabilities from each subtype
smr_liu <- do.call(
  rbind,
  bplapply(subtype_dirs, function(.x) wrapper1(.x))
)
print("Finished processing Liu2021 data.")

tp_liu <- vmrseq:::.estimTransitProbsFromSummary(
  smr_units = smr_liu,
  max_dist_bp = max_dist_bp,
  buffer_bp = buffer_bp,
  degree = degree,
  span = span
)
plotTransitProbs(tp_liu)
ggsave(paste0("plots/estim_transitProb_subtypes/estimated_transit_probs_liu_lag",
              lags[1],"_",lags[length(lags)],".png"), 
       width = 8, height = 6)

##### Sample cells in Luo2017 mice&human data and form into list of data.frames ====

wrapper2 <- function(subtype_dir){
  
  st <- loadHDF5SummarizedExperiment(subtype_dir) %>% 
    granges() %>% 
    as.data.table() %>%
    mutate(state = round(cell_MF)) %>%
    group_by(seqnames) %>%
    mutate(lag_state = lag(state, 1))
  
  # decide the states for sites with MF = 0.5 based on its previous site
  st$state[which(st$cell_MF==0.5)] <- st$lag_state[which(st$cell_MF==0.5)]
  
  smr <- vmrseq:::.computeProb1Unit(
    df = st %>% select(seqnames, start, state),
    max_dist_bp = max_dist_bp,
    buffer_bp = buffer_bp,
    lags = lags
  )
  
  return(smr)
}

# Luo2017 mice data

dir <- ("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/data/processed/processed_luo2017_mice/")
subtype_dirs <- paste0(dir, grep(pattern = "qced", x = list.files(dir), value = T))  # QCed: across-cell coverage >= 5

smr_luom <- do.call(
  rbind,
  bplapply(subtype_dirs, function(.x) wrapper2(.x))
)
print("Finished processing Luo2017 mice data.")

tp_luom <- vmrseq:::.estimTransitProbsFromSummary(
  smr_units = smr_luom,
  max_dist_bp = max_dist_bp,
  buffer_bp = buffer_bp,
  degree = degree,
  span = span
)
plotTransitProbs(tp_luom)
ggsave(paste0("plots/estim_transitProb_subtypes/estimated_transit_probs_luom_lag",
              lags[1],"_",lags[length(lags)],".png"), 
       width = 8, height = 6)


# Luo2017 human data

dir <- ("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/data/processed/processed_luo2017_human/")
subtype_dirs <- paste0(dir, grep(pattern = "qced", x = list.files(dir), value = T))  # QCed: across-cell coverage >= 5

smr_luoh <- do.call(
  rbind,
  bplapply(subtype_dirs, function(.x) wrapper2(.x))
)
print("Finished processing Luo2017 human data.")

tp_luoh <- vmrseq:::.estimTransitProbsFromSummary(
  smr_units = smr_luoh,
  max_dist_bp = max_dist_bp,
  buffer_bp = buffer_bp,
  degree = degree,
  span = span
)
plotTransitProbs(tp_luoh)
ggsave(paste0("plots/estim_transitProb_subtypes/estimated_transit_probs_luoh_lag",
              lags[1],"_",lags[length(lags)],".png"), 
       width = 8, height = 6)

#########################################################################
##### Estimate transition probability from all cells from 3 datasets ####
#########################################################################

smr <- rbind(smr_liu, smr_luom, smr_luoh)
saveRDS(smr, "code/estim_transitProb/summary.rds")

tp0 <- vmrseq:::.estimTransitProbsFromSummary(
  smr_units = smr,
  max_dist_bp = max_dist_bp,
  buffer_bp = buffer_bp,
  degree = degree,
  span = span
)

# plot loess-fitted transition probs
plotTransitProbs(tp0)
ggsave(paste0("plots/estim_transitProb_subtypes/estimated_transit_probs_lag",
              lags[1],"_",lags[length(lags)],".png"), 
       width = 8, height = 6)

# Save tp0 object
saveRDS(tp0, "code/estim_transitProb/tp0.rds")