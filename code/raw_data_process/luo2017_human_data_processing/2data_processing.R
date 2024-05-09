source("code/SETPATHS.R")
source("code/luo2017_human_data_processing/helper_functions.R")

#### NOTICE: Cells are selected according to Table S2 reported by scMET paper (partial cells removed compared to Luo2017 Table S1)
#### scMET ref: Kapourani et al. scMET: Bayesian modeling of DNA methylation heterogeneity at single-cell resolution. Genome Biol 22, 114 (2021). https://doi.org/10.1186/s13059-021-02329-8)
#### Luo2017 ref: Luo et al. Single-cell methylomes identify neuronal subtypes and regulatory elements in mammalian cortex. Science. 2017 Aug 11;357(6351):600-604. doi: 10.1126/science.aan3351. PMID: 28798132; PMCID: PMC5570439.
metadata <- fread("data/metadata/metadata_luo2017/NIHMS893063-supplement-Table_S2_csv.csv", skip = 1)[, .(Sample, `Neuron type`)]
metadata[, .(.N), by = `Neuron type`] %>% arrange(desc(N)) %>% filter(N >= 100)
#    Neuron type   N
# 1:       hL2/3 873
# 2:       hPv-1 175
# 3:       hNdnf 173
# 4:      hSst-2 163
# 5:       hL5-4 162
# 6:       hDL-1 144
# 7:      hVip-2 128
# 8:         hL4 109
# 9:       hL6-2 107

# wrapper(subtype = "hL2/3", mc.cores = 16) 
# wrapper(subtype = "hPv-1", mc.cores = 16)
# wrapper(subtype = "hNdnf", mc.cores = 16)
# wrapper(subtype = "hSst-2", mc.cores = 16)
# wrapper(subtype = "hL5-4", mc.cores = 16)
# wrapper(subtype = "hDL-1", mc.cores = 16)
# wrapper(subtype = "hVip-2", mc.cores = 16)
# wrapper(subtype = "hL4", mc.cores = 16)
wrapper(subtype = "hL6-2", mc.cores = 8)




