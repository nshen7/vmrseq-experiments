## Objective of this script: 
## To estimate (default) transition probability distribution from subsampled cells in (unQC'ed) subtypes used for training emission beta-mixture priors. 

library(tidyverse)
setwd(here::here())
source("code/package_functions/transitProb.R")

n_cells <- 50 # take 50 cells from each subtype for estimation

################################################################
# ==== Compute cell summaries in Liu2021 data ====
################################################################

wrapper1 <- function(info, n_cells, max_dist_bp, buffer_bp, seed = 2022){
  dir <- "data/processed/processed_liu2021/"
  file_list <- list.files(dir)
  file_name <- grep(pattern = paste0("sample", info$Sample,
                                     "_", info$GEO_accession, 
                                     "_subtype_", sub(" ", "_", info$SubType),
                                     ".*cells$"),
                    x = file_list, value = T)
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))
  
  set.seed(seed)
  index <- sample(ncol(cells.se), n_cells)
  smr_cells <- do.call(rbind, 
                       map(index, ~ .computeProb1Cell(df = data.frame(chr = seqnames(cells.se),
                                                                      pos = start(cells.se),
                                                                      MF = round(assays(cells.se[, .x])$M / assays(cells.se[, .x])$Cov)) %>% na.omit,
                                                      max_dist_bp = max_dist_bp, buffer_bp = buffer_bp)
                           )
                       )
  return(smr_cells)
}

## information of all subtypes used for training
metadata <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv")
cell_counts <- metadata[, .(.N), by = .(Sample, CellClass, SubType, GEO_accession)]

cell_counts_exc <- cell_counts %>% # selected subtypes in excitatory cells for training
  filter(CellClass=="Exc") %>%
  arrange(desc(N)) %>%
  filter(!is.na(SubType), N >= 100) %>%
  dplyr::slice(c(5L, which(!duplicated(SubType)&!duplicated(Sample)))) %>%
  dplyr::slice(-2L)

cell_counts_inh <- cell_counts %>%
  filter(CellClass=="Inh") %>%
  arrange(desc(N)) %>%
  filter(!is.na(SubType), N >= 100) %>%
  dplyr::slice(which(!duplicated(SubType)&!duplicated(Sample)))

info <- rbind(cell_counts_exc,cell_counts_inh)

## Compute empirical probabilities from individual cells
smr_liu2021 <- do.call(rbind, parallel::mclapply(1:nrow(info), 
                                                 function(.x) wrapper1(info[.x,], n_cells = n_cells, max_dist_bp = 2000, buffer_bp = 3000),
                                                 mc.cores = 16))
fwrite(smr_liu2021, paste0("data/interim/estim_transitProb/transitProb_",nrow(info),"subtypes_",nrow(info)*n_cells,"cells_summary_Liu2021.csv"))



################################################################
# ==== Compute cell summaries in Luo2017 data ====
################################################################

wrapper2 <- function(info, n_cells, max_dist_bp, buffer_bp, seed = 2022){
  dir <- "data/processed/processed_luo2017_mice/"
  file_list <- list.files(dir)
  file_name <- grep(pattern = paste0("subtype_", sub("/", "", info$SubType), ".*cells$"),
                    x = file_list, value = T)
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))
  
  set.seed(seed)
  index <- sample(ncol(cells.se), n_cells)
  smr_cells <- do.call(rbind, 
                       map(index, ~ .computeProb1Cell(df = data.frame(chr = seqnames(cells.se),
                                                                      pos = start(cells.se),
                                                                      MF = round(assays(cells.se[, .x])$M / assays(cells.se[, .x])$Cov)) %>% na.omit,
                                                      max_dist_bp = max_dist_bp, buffer_bp = buffer_bp)
                       )
  )
  return(smr_cells)
}

## information of all subtypes used for training
metadata <- fread("../../DXM_extend_chr1/data/metadata/sample_info_processed.csv")
info <- metadata[, .(.N), by = .(Neuron_type1, Neuron_type3)] %>% 
  arrange(desc(N)) %>%
  filter(N >= 100) %>%
  dplyr::rename(CellClass = Neuron_type1, SubType = Neuron_type3) %>%
  mutate(CellClass = recode(CellClass, "Excitatory" = "Exc", "Inhibitory" = "Inh"))

## compute density and fit parameters
smr_luo2017 <- do.call(rbind, parallel::mclapply(1:nrow(info), 
                                                 function(.x) wrapper2(info[.x, ], n_cells = n_cells, max_dist_bp = 2000, buffer_bp = 3000),
                                                 mc.cores = 16))
fwrite(smr_luo2017, paste0("data/interim/estim_transitProb/transitProb_",nrow(info),"subtypes_",nrow(info)*n_cells,"cells_summary_Luo2017.csv"))



################################################################
# ==== Estimate transition probability from 27 subtypes ====
################################################################

smr_cells <- rbind(fread("data/interim/estim_transitProb/transitProb_18subtypes_900cells_summary_Liu2021.csv"), 
                   fread("data/interim/estim_transitProb/transitProb_9subtypes_450cells_summary_Luo2017.csv"))

tp <- .estimTransitProbsFromSummary(smr_cells, max_dist_bp = 2000, buffer_bp = 3000, cores = 16, degree = 2, span = 0.02)
write_rds(tp, here::here("code/package_functions/transitProbs_27subtypes_1350cells_Luo2017&Liu2021.rds"))

plotTransitProbs(tp, fitted = T, linewidth = 0.75) 
ggsave("code/package_functions/transitProbs_27subtypes_270cells_Luo2017&Liu2021.png", width = 7, height = 5)


