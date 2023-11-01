source("code/SETPATHS.R")
write_dir <- here("manuscript_related", "manuscript_tables")

md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv") %>%
  select(-c(Neuron_type, Neuron_type2, Neuron_type3)) %>%
  rename(BroadClass = Neuron_type1, Subtype = Neuron.type)
fwrite(md, here(write_dir, 'supp_table_sample_metadata_luo.csv'))
