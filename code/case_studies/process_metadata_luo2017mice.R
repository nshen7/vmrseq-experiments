library(data.table)
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") 

# Original metadata table from Luo2017
md_0 <- read.csv("data/metadata/metadata_luo2017/luo_supp_tables/NIHMS893063-supplement-Table_S1_csv.csv", skip = 1) %>%
  mutate(Sample = gsub("_indexed", "", Sample))

# Processed metadata table from scMET paper
md_1 <- as.vector(t(read.csv("data/metadata/metadata_luo2017/scmet_supp_tables/Table_S2_Luo2017_sample_metadata.csv",
                                        stringsAsFactors = F, header=FALSE)))
sample_info <- do.call(rbind.data.frame, sapply(md_1[-1], function(x){strsplit(x,split="\t")}))
colnames(sample_info) <- strsplit(md_1[1], "\t")[[1]]

# Add original Neuron.type info to md_1
sample_info <- sample_info %>% 
  left_join(md_0[c('Sample','Neuron.type')], by = c('sample' = 'Sample')) 

write.csv(sample_info, "data/metadata/metadata_luo2017/sample_info_processed.csv", quote = F, row.names = F)
