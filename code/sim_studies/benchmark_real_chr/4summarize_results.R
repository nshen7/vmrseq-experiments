.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
devtools::load_all("../vmrseq-package/vmrseq/")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(GenomicRanges))

N <- 200
NP <- 4
NV <- 5000
seed <- 2022

# ==== vmrseq ====

vmrseq_eval <- function(res_vmrseq, penalty) {
  res_gr <- res_vmrseq$gr_qced
  index <- which(res_gr$loglik_diff < penalty) # indices of sites failed to pass penalty threshold
  res_gr$vmr_index[index] <- NA
  res_gr$vmr_num_cpg[index] <- NA
  power_site <- with(res_gr, sum(is_vml&!is.na(vmr_index))/sum(is_vml))
  fdr_site <- with(res_gr, sum(!is_vml&!is.na(vmr_index))/sum(!is.na(vmr_index)))
  # cat("Site-level power =", power_site, "; ")
  # cat("Site-level FDR =", fdr_site, "\n")
  
  true <- res_gr %>% data.frame %>%
    group_by(cr_index) %>%
    mutate(cr_n_cpg = sum(!is.na(cr_index))) %>%
    ungroup() %>%
    filter(is_vml) %>%
    group_by(vmr_name) %>%
    summarise(n_cpg = n(),
              n_dectected = sum(!is.na(vmr_index)),
              n_incr = sum(!is.na(cr_index)),
              cr_index = max(cr_index, na.rm = T),
              cr_n_cpg = max(cr_n_cpg),
              pi1 = max(pi1, na.rm = T),
              med_total = median(total),
              loglik_diff = max(loglik_diff, na.rm = T)) %>%
    filter(n_cpg >= 5)
  power_region <- sum(true$n_dectected > 0) / nrow(true)
  # cat("Region-level power =", power, "; ")
  detected <- res_gr %>% data.frame %>%
    filter(!is.na(vmr_index)) %>%
    group_by(cr_index) %>%
    mutate(cr_n_cpg = n()) %>%
    ungroup() %>%
    group_by(vmr_index) %>%
    summarise(n_cpg = n(),
              n_true = sum(!is.na(vmr_name)),
              cr_index = max(cr_index),
              cr_n_cpg = max(cr_n_cpg),
              pi1 = max(pi1, na.rm = T),
              optim_pi = max(pi),
              mean_MF = mean(meth/total),
              loglik_diff = max(loglik_diff))
  fdr_region <- 1-sum(detected$n_true > 0) / nrow(detected)
  # cat("Region-level FDR =", fdr, "\n")
  return(list(power_site = power_site, fdr_site = fdr_site, 
              power_region = power_region , fdr_region  = fdr_region , 
              true = true, detected = detected))
}

res_vmrseq <- readRDS(
  paste0(
    "data/interim/sim_studies/benchmark_real_chr/vmrseq/output/pseudoChr_",
    subtype, "_", chromosome, "_", 
    N, "cells_", NP, "subpops_", 
    NV, "VMRs_seed", seed, "_vmrseqOutput.rds"
  )
)

smr_vmrseq <- data.frame(method = "vmrseq",
                         penalty = c(0:10, seq(15,60,5)),
                         power_site = NA,
                         fdr_site = NA,
                         power_region = NA,
                         fdr_region = NA)

for (i in 1:nrow(smr_vmrseq)) {
  pen <- smr_vmrseq$penalty[i]
  smr_temp <- vmrseq_eval(res_vmrseq, penalty = pen)
  if (pen==0) {View(smr_temp$true); View(smr_temp$detected)}
  smr_vmrseq$penalty[i] <- pen
  smr_vmrseq$power_site[i] <- smr_temp$power_site
  smr_vmrseq$fdr_site[i] <- smr_temp$fdr_site
  smr_vmrseq$power_region[i] <- smr_temp$power_region
  smr_vmrseq$fdr_region[i] <- smr_temp$fdr_region
  cat(i, " ")
}

# with(smr_vmrseq, plot(fdr_site, power_site, type = "l", xlim = c(0,1), ylim = c(0,1)))
# with(smr_vmrseq, plot(fdr_region, power_region, type = "l", xlim = c(0,1), ylim = c(0,1)))

# === scbs ====
gr <- loadHDF5SummarizedExperiment(
  paste0(
    "data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
    subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
    NV, "VMRs_seed", seed
  )
) %>% granges()

smr_scbs <- data.frame(method = "scbs",
                       vt = c(seq(.001, .009, .001), seq(.01, .09, .01), seq(0.1,0.3,0.1)),
                       # vt = c(seq(.001, .009, .001), seq(.01, .08, .01)),
                       power_site = NA,
                       fdr_site = NA,
                       power_region = NA,
                       fdr_region = NA)

for (i in 1:nrow(smr_scbs)) {
  vt <- smr_scbs$vt[i]
  res_scbs <- fread(
    paste0("data/interim/sim_studies/benchmark_real_chr/scbs/output/pseudoChr_IT-L23_Cux1_chr1_", 
           N, "cells_", NP, "subpops_", 
           NV, "VMRs_seed", seed, "_", vt, "vt.bed")
  )
  colnames(res_scbs) <- c("chr", "start", "end", "meth_var")
  
  res_scbs.gr <- with(
    res_scbs, 
    GRanges(seqnames = chr,
            ranges = IRanges(start = start, end = end))
  )
  
  hits <- findOverlaps(gr, res_scbs.gr)
  gr$vmr_index <- NA # to store detected VMRs by scbs
  gr$vmr_index[queryHits(hits)] <- subjectHits(hits) 
  smr_scbs$power_site[i] <- with(gr, sum(is_vml&!is.na(vmr_index)) / sum(is_vml))
  smr_scbs$fdr_site[i] <- with(gr, sum(!is_vml&!is.na(vmr_index)) / sum(!is.na(vmr_index)))
  
  true_olap <- gr %>% as.data.frame %>% group_by(vmr_name) %>%  
    summarise(n_olap = sum(!is.na(vmr_index))) %>% select(n_olap) %>% unlist
  smr_scbs$power_region[i] <- sum(true_olap>0)/length(true_olap)
  
  detected_olap <- gr %>% as.data.frame %>% group_by(vmr_index) %>%  
    summarise(n_olap = sum(!is.na(vmr_name))) %>% select(n_olap) %>% unlist
  smr_scbs$fdr_region[i] <- sum(detected_olap==0)/length(detected_olap)
    
  
  cat(i, " ")
}

# with(smr_scbs, lines(fdr_site, power_site, col = "red"))
# with(smr_scbs, lines(fdr_region, power_region, col = "red"))

# ==== scmet ====
bp_size <- 20000
gr <- loadHDF5SummarizedExperiment(
  paste0(
    "data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
    subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
    NV, "VMRs_seed", seed
  )
) %>% granges()
wds.gr <- readRDS(paste0(
  "data/interim/sim_studies/benchmark_real_chr/scmet/input/features_", 
  subtype, "_", chromosome, "_", 
  bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
  NV, "VMRs_seed", seed, ".rds"
))

smr_scmet <- data.frame(method = paste0("scMET_", bp_size/1000, "kb"),
                       efdr = c(0.01, 0.02, 0.05, seq(0.1,0.9,0.1), 0.99),
                       power_site = NA,
                       fdr_site = NA,
                       power_region = NA,
                       fdr_region = NA)
for (i in 1:nrow(smr_scmet)) {
  
  efdr <- smr_scmet$efdr[i]
  res_scmet <- fread(paste0(
    "data/interim/sim_studies/benchmark_real_chr/scmet/output/pseudoChr_",
    subtype, "_", chromosome, "_", 
    bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
    NV, "VMRs_seed", seed, "_", efdr, "efdr_scmetOutput.csv"
  )) %>% mutate(window = as.integer(gsub(".*_(.*)$", "\\1", feature_name))) %>% filter(is_variable)
  
  hvf.gr <- wds.gr[res_scmet$window]
  
  hits <- findOverlaps(gr, hvf.gr)
  gr$vmr_index <- NA # to store detected VMRs by scbs
  gr$vmr_index[queryHits(hits)] <- subjectHits(hits) 
  smr_scmet$power_site[i] <- with(gr, sum(is_vml&!is.na(vmr_index)) / sum(is_vml))
  smr_scmet$fdr_site[i] <- with(gr, sum(!is_vml&!is.na(vmr_index)) / sum(!is.na(vmr_index)))
  
  true_olap <- gr %>% as.data.frame %>% group_by(vmr_name) %>%  
    summarise(n_olap = sum(!is.na(vmr_index))) %>% select(n_olap) %>% unlist
  smr_scmet$power_region[i] <- sum(true_olap>0)/length(true_olap)
  
  detected_olap <- gr %>% as.data.frame %>% group_by(vmr_index) %>%  
    summarise(n_olap = sum(!is.na(vmr_name))) %>% select(n_olap) %>% unlist
  smr_scmet$fdr_region[i] <- sum(detected_olap==0)/length(detected_olap)
  
  
}

# ==== plots ====
smr <- rbind(smr_vmrseq[-2], smr_scbs[-2], smr_scmet[-2])
# smr <- rbind(smr_vmrseq[-2], smr_scbs[-2])

# site-level fdr and power
smr %>%  
  ggplot(aes(fdr_site, power_site, color = method)) + 
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dashed") +
  geom_path() +
  xlab("Site-level FDR") + ylab("Site-level power") +
  ggtitle("Modified real chromosome (200 cells, 4 subpops)") +
  xlim(0,1) + ylim(0,1)
ggsave(paste0(
  "plots/sim_studies/benchmark_real_chr/comparison/pseudoChr_fdr&power_siteLevel_", 
  N , "cells_", NP, "subpops_seed", seed, ".png"
), width = 6, height = 5)

# region-level fdr and power
smr %>%  
  ggplot(aes(fdr_region, power_region, color = method)) + 
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dashed") +
  geom_path() +
  xlab("Region-level FDR") + ylab("Region-level power") +
  ggtitle("Modified real chromosome (200 cells, 4 subpops)") +
  xlim(0,1) + ylim(0,1)
ggsave(paste0(
  "plots/sim_studies/benchmark_real_chr/comparison/pseudoChr_fdr&power_regionLevel_", 
  N , "cells_", NP, "subpops_seed", seed, ".png"
), width = 6, height = 5)
