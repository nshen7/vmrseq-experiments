source("code/SETPATHS.R")
write_dir <- here("manuscript_related", "manuscript_tables")

md <- fread(here('data', 'metadata', 'argelaguet2019', 'argelaguet2019_full_met&rna_sample_metadata_processed.csv'))
fwrite(md, here(write_dir, 'supp_table_sample_metadata_argelaguet.csv'))

## Source of 'supp_table_argelaguet_gene_vmr_corr.csv' is in:
##  '/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/manuscript_related/manuscript_figures/figs_argelaguet_study.R'