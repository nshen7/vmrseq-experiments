source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")

write_dir <- here("manuscript_related", "manuscript_tables")

# ---- read in all subtypes ----

# According to code/raw_data_process/liu2021_data_processing/2data_processing_subtypeAcrossSample.R
md_liu <- fread("data/metadata/metadata_liu2021/Liu2021_cell_full_metadata_processed.csv") %>%
  filter(!is.na(GEO_accession) & !is.na(SubType) & !is.na(FilePath))
st_liu <- md_liu[, .(.N), by = .(CellClass, SubType)] %>% 
  filter(!grepl("Outlier", SubType)) %>%
  mutate(dataset = 'Liu et al. 2021')

# According to code/(archived)estim_emiBetaPrior/1data_subsampling.R
md_luo_m <- fread("../../DXM_extend_chr1/data/metadata/sample_info_processed.csv")
st_luo_m <- md_luo_m[, .(.N), by = .(Neuron_type1, Neuron_type3)] %>% 
  dplyr::rename(CellClass = Neuron_type1, SubType = Neuron_type3) %>%
  mutate(CellClass = recode(CellClass, "Excitatory" = "Exc", "Inhibitory" = "Inh")) %>%
  mutate(dataset = 'Luo et al. 2017 (Mouse)')

# According to code/(archived)estim_emiBetaPrior/1data_subsampling.R
md_luo_h <- fread("data/metadata/metadata_luo2017/luo_supp_tables/NIHMS893063-supplement-Table_S2_csv.csv", skip = 1)[, .(Sample, `Neuron type`)]
st_luo_h <- md_luo_h[, .(.N), by = .(`Neuron type`)] %>%
  dplyr::rename(SubType = `Neuron type`)  %>%
  mutate(CellClass = 'NA', dataset = 'Luo et al. 2017 (Human)')

# ---- add in columns of indicator for used in training transition probs ----

# According to code/raw_data_process/liu2021_data_processing/2data_processing_subtypeAcrossSample.R
st_sel_liu <- fread("data/metadata/metadata_liu2021/Liu2021_subtypes_cellCount_binnedBy100CellCount_selectedForParamTraining.csv")
st_liu <- st_liu %>% mutate(used_for_train_transition_prob = (SubType %in% st_sel_liu$SubType))

# According to code/raw_data_process/luo2017_mice_data_processing/2data_processing.R
st_luo_m <- st_luo_m %>% mutate(used_for_train_transition_prob = (N >= 100))

# According to code/raw_data_process/luo2017_human_data_processing/2data_processing.R
st_luo_h <- st_luo_h %>% mutate(used_for_train_transition_prob = (N >= 100))


# ---- add in columns of indicator for used in training beta priors in emission probs ----

# According to /scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-package/vmrseq/data-raw/PARAMS.R
sub_liu2021 <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_liu2021_clusterMeth.csv")
sub_luo2017m <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017mice_clusterMeth.csv")
sub_luo2017h <- fread("data/interim/estim_emiBetaPrior_ZIBBregression2/emiBetaPrior_subtype_subsample_luo2017human_clusterMeth.csv")

st_liu <- st_liu %>% mutate(used_for_train_beta_prior = (SubType %in% sub_liu2021$SubType))
st_luo_m <- st_luo_m %>% mutate(used_for_train_beta_prior = (SubType %in% sub_luo2017m$SubType))
st_luo_h <- st_luo_h %>% mutate(used_for_train_beta_prior = (SubType %in% sub_luo2017h$SubType))


# ---- write out ----

st_smr <- rbind(st_liu, st_luo_m, st_luo_h) %>% 
  rename('cell_class' = 'CellClass',
         'subtype' = 'SubType',
         'cell_count' = 'N')
fwrite(st_smr, here(write_dir, "supp_table_subtypes_in_param_training.csv"))

