suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(data.table))

setwd(here::here())

## number of sites sampled from each subtype
sample_n <- 50000
set.seed(2022)

################################################################
# ==== Subsample subtypes in Liu2021 data ====
################################################################
wrapper1 <- function(info_row, cores = 16){
  dir <- "data/processed/processed_liu2021/"
  file_list <- list.files(dir)
  file_name <- grep(pattern = paste0("sample", info_row$Sample,
                                     "_", info_row$GEO_accession,
                                     "_subtype_", sub(" ", "_", info_row$SubType),
                                     ".*_qced$"),
                    x = file_list, value = T)
  print(file_name)
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))
  
  index <- sample(1:nrow(cells.se), sample_n)
  sub_dat <- data.frame(
    DataSource = "Liu2021",
    CellClass = info_row$CellClass,
    SubType = info_row$SubType,
    N_cell = ncol(cells.se),
    med_cov = median(values(cells.se)$cell_cov),
    cell_cov = values(cells.se)$cell_cov[index],
    cell_meth = values(cells.se)$cell_meth[index],
    cell_MF = values(cells.se)$cell_MF[index]
  )
  
  return(sub_dat)
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

## compute density and fit parameters
sub_liu2021 <- do.call(rbind, map(1:nrow(info), function(i) wrapper1(info[i])))
fwrite(sub_liu2021, "data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_liu2021.csv")


################################################################
# ==== Subsample subtypes in Luo2017 data ====
################################################################

wrapper2 <- function(info_row, cores = 16){
  dir <- "data/processed/processed_luo2017_mice/"
  file_list <- list.files(dir)
  file_name <- grep(pattern = paste0("subtype_", sub("/", "", info_row$SubType), ".*_qced$"),
                    x = file_list, value = T)
  print(file_name)
  cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))
  
  index <- sample(1:nrow(cells.se), sample_n)
  sub_dat <- data.frame(
    DataSource = "Luo2017",
    CellClass = info_row$CellClass,
    SubType = info_row$SubType,
    N_cell = ncol(cells.se),
    med_cov = median(values(cells.se)$cell_cov),
    cell_cov = values(cells.se)$cell_cov[index],
    cell_meth = values(cells.se)$cell_meth[index],
    cell_MF = values(cells.se)$cell_MF[index]
  )
  
  return(sub_dat)
}

## information of all subtypes used for training
metadata <- fread("../../DXM_extend_chr1/data/metadata/sample_info_processed.csv")
info <- metadata[, .(.N), by = .(Neuron_type1, Neuron_type3)] %>%
  arrange(desc(N)) %>%
  filter(N >= 100) %>%
  dplyr::rename(CellClass = Neuron_type1, SubType = Neuron_type3) %>%
  mutate(CellClass = recode(CellClass, "Excitatory" = "Exc", "Inhibitory" = "Inh"))

## compute density and fit parameters
sub_luo2017 <- do.call(rbind, map(1:nrow(info), function(i) wrapper2(info[i])))
fwrite(sub_luo2017, "data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_luo2017.csv")


#############################################
#### Processed into meth and unmeth data ####
#############################################

sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_liu2021.csv")
sub_luo2017 <- fread("data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_luo2017.csv")
sub_unmeth <- rbind(sub_luo2017, sub_liu2021) %>% 
  as_tibble() %>%
  dplyr::filter(cell_MF < 0.5) %>%
  mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth)))
saveRDS(sub_unmeth, "data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_unmethClust.rds")
sub_meth <- rbind(sub_luo2017, sub_liu2021) %>% 
  as_tibble() %>%
  dplyr::filter(cell_MF > 0.5) %>%
  mutate(y = as.matrix(data.frame(cell_meth, cell_cov - cell_meth)))
saveRDS(sub_meth, "data/interim/estim_emiBetaPrior_ZIBBregression/emiBetaPrior_subtype_subsample_methClust.rds")
