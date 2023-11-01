source("code/SETPATHS.R")
devtools::load_all('../vmrseq-package/vmrseq/')

write_dir <- here("manuscript_related", "manuscript_tables")

## Default transition probabilities
tp0 <- vmrseq:::tp0@transit_probs %>%
  cbind(distance_bp = 1:nrow(.), .)
fwrite(tp0, here(write_dir, 'supp_table_default_transit_probs.csv'))

## Default beta priors in emission probabilities
params_m <- vmrseq:::params_m
fwrite(params_m, here(write_dir, 'supp_table_default_beta_priors_m.csv'))
params_u <- vmrseq:::params_u %>%
  rename(median_coverage = med_cov)
fwrite(params_u, here(write_dir, 'supp_table_default_beta_priors_u.csv'))
