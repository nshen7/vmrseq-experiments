setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
source("code/luo2017_mice_data_processing/helper_functions.R")

#### NOTICE: Cells are selected according to Table S2 reported by scMET paper (partial cells removed compared to Luo2017 Table S1)
#### scMET ref: Kapourani et al. scMET: Bayesian modeling of DNA methylation heterogeneity at single-cell resolution. Genome Biol 22, 114 (2021). https://doi.org/10.1186/s13059-021-02329-8)
#### Luo2017 ref: Luo et al. Single-cell methylomes identify neuronal subtypes and regulatory elements in mammalian cortex. Science. 2017 Aug 11;357(6351):600-604. doi: 10.1126/science.aan3351. PMID: 28798132; PMCID: PMC5570439.

metadata <- fread("../../DXM_extend_chr1/data/metadata/sample_info_processed.csv")
cell_counts <- metadata[, .(.N), by = .(Neuron_type1, Neuron_type3)] %>% arrange(desc(N)) 
#     Neuron_type1 Neuron_type3   N
#  1:   Excitatory        mL6-2 686 (processed)
#  2:   Excitatory        mL2/3 649 (processed)
#  3:   Excitatory          mL4 370 (processed)
#  4:   Excitatory        mL5-1 284 (processed)
#  5:   Excitatory        mDL-2 272 (processed)
#  6:   Inhibitory          mPv 136 (processed)
#  7:   Excitatory        mL5-2 128
#  8:   Inhibitory      mSst-12 123 (processed)
#  9:   Excitatory        mDL-1 122
# 10:   Inhibitory     mNdnf-12  96
# 11:   Excitatory        mL6-1  80
# 12:   Inhibitory         mVip  69
# 13:   Excitatory        mDL-3  54



# ==== Exitatory ====
# ### Neuron_type3 = mDL-2; n.cells = 272 finished
# Sys.time()
# wrapper(subtype = "mDL-2", mc.cores = 16)
# Sys.time() 

# ### Neuron_type3 = mL4; n.cells = 370 finished
# Sys.time()
# wrapper(subtype = "mL4", mc.cores = 16) 
# Sys.time()

### Neuron_type3 = mL2/3; n.cells = 649 finished
# Sys.time()
# wrapper(subtype = "mL2/3", mc.cores = 16)
# Sys.time()

# ### Neuron_type3 = mL5-1; n.cells = 284 
# Sys.time()
# wrapper(subtype = "mL5-1", mc.cores = 16)
# Sys.time()

# ### Neuron_type3 = mL5-2; n.cells = 128
# Sys.time()
# wrapper(subtype = "mL5-2", mc.cores = 16)
# Sys.time()

### Neuron_type3 = mDL-1; n.cells = 122
Sys.time()
wrapper(subtype = "mDL-1", mc.cores = 16)
Sys.time()

# ==== Inhibitory ====
# ### Neuron_type3 = mPv; n.cells = 136
# Sys.time()
# wrapper(subtype = "mPv", mc.cores = 16)
# Sys.time()

# ### Neuron_type3 = mSst-12; n.cells = 123
# Sys.time()
# wrapper(subtype = "mSst-12", mc.cores = 16)
# Sys.time()




