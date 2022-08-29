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
seed <- 2022
chromosome <- "chr1"

# ==== vmrseq ====

.vmrseqEval <- function(res_vseq, penalty) {
  res_gr <- res_vseq$gr_qced
  index <- which(res_gr$loglik_diff < penalty) # indices of sites failed to pass penalty threshold
  res_gr$vmr_index[index] <- NA
  res_gr$vmr_num_cpg[index] <- NA
  power_site <- with(res_gr, sum(is_vml&!is.na(vmr_index))/sum(is_vml))
  fdr_site <- with(res_gr, sum(!is_vml&!is.na(vmr_index))/sum(!is.na(vmr_index)))
  # cat("Site-level power =", power_site, "; ")
  # cat("Site-level FDR =", fdr_site, "\n")
  
  true <- res_gr %>% data.frame %>%
    filter(is_vml) %>%
    group_by(vmr_name) %>%
    summarise(n_cpg = n(),
              n_dectected = sum(!is.na(vmr_index)),
              n_incr = sum(!is.na(cr_index)),
              cr_index = max(cr_index, na.rm = T),
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
    summarise(cr_index = max(cr_index),
              n_cpg = n(),
              n_true = sum(!is.na(vmr_name)),
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

.vmrseqCrEval <- function(res_vseq, cutoff) {
  res_gr <- res_vseq$gr_qced
  index <- which(res_gr$var < cutoff) # indices of sites failed to pass penalty threshold
  res_gr$cr_index[index] <- NA
  power_site <- with(res_gr, sum(is_vml&!is.na(cr_index))/sum(is_vml))
  fdr_site <- with(res_gr, sum(!is_vml&!is.na(cr_index))/sum(!is.na(cr_index)))
  # cat("Site-level power =", power_site, "; ")
  # cat("Site-level FDR =", fdr_site, "\n")
  
  true <- res_gr %>% data.frame %>%
    filter(is_vml) %>%
    group_by(vmr_name) %>%
    summarise(n_cpg = n(),
              n_dectected = sum(!is.na(cr_index)),
              n_incr = sum(!is.na(cr_index)),
              cr_index = max(cr_index, na.rm = T),
              pi1 = max(pi1, na.rm = T),
              med_total = median(total),
              loglik_diff = max(loglik_diff, na.rm = T)) %>%
    filter(n_cpg >= 5)
  power_region <- sum(true$n_dectected > 0) / nrow(true)
  # cat("Region-level power =", power, "; ")
  detected <- res_gr %>% data.frame %>%
    filter(!is.na(cr_index)) %>%
    group_by(cr_index) %>%
    mutate(cr_n_cpg = n()) %>%
    ungroup() %>%
    group_by(cr_index) %>%
    summarise(cr_index = max(cr_index),
              n_cpg = n(),
              n_true = sum(!is.na(vmr_name)),
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

summarizeResVseq <- function(NV) {
  res_vseq <- readRDS(
    paste0(
      "data/interim/sim_studies/benchmark_real_chr/vmrseq/output/pseudoChr_",
      subtype, "_", chromosome, "_", 
      N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed, "_vmrseqOutput.rds"
    )
  )
  
  smr_vseq <- data.frame(method = "vmrseq",
                         NV = NV,
                         penalty = c(-2:10, seq(15,60,5)),
                         power_site = NA,
                         fdr_site = NA,
                         power_region = NA,
                         fdr_region = NA)
  
  for (i in 1:nrow(smr_vseq)) {
    pen <- smr_vseq$penalty[i]
    smr <- .vmrseqEval(res_vseq, penalty = pen)
    # if (i==1) {View(smr$true); View(smr$detected)}
    smr_vseq$penalty[i] <- pen
    smr_vseq$power_site[i] <- smr$power_site
    smr_vseq$fdr_site[i] <- smr$fdr_site
    smr_vseq$power_region[i] <- smr$power_region
    smr_vseq$fdr_region[i] <- smr$fdr_region
    cat(i, " ")
  }
  return(smr_vseq)
}

summarizeResVseqCr <- function(NV) {
  res_vseq <- readRDS(
    paste0(
      "data/interim/sim_studies/benchmark_real_chr/vmrseq/output/pseudoChr_",
      subtype, "_", chromosome, "_", 
      N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed, "_vmrseqOutput.rds"
    )
  )
  smr_vseq_cr <- data.frame(method = "vmrseq_CR",
                            NV = NV,
                            cutoff = seq(0.05,0.2,0.01),
                            power_site = NA,
                            fdr_site = NA,
                            power_region = NA,
                            fdr_region = NA)
  
  for (i in 1:nrow(smr_vseq_cr)) {
    cutoff <- smr_vseq_cr$cutoff[i]
    smr <- .vmrseqCrEval(res_vseq, cutoff = cutoff)
    smr_vseq_cr$cutoff[i] <- cutoff
    smr_vseq_cr$power_site[i] <- smr$power_site
    smr_vseq_cr$fdr_site[i] <- smr$fdr_site
    smr_vseq_cr$power_region[i] <- smr$power_region
    smr_vseq_cr$fdr_region[i] <- smr$fdr_region
    cat(i, " ")
  }
  return(smr_vseq_cr)
}

# === scbs ====

summarizeResScbs <- function(NV) {
  gr <- loadHDF5SummarizedExperiment(
    paste0(
      "data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
      subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed
    )
  ) %>% granges()
  
  smr_scbs <- data.frame(method = "scbs",
                         NV = NV,
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
  return(smr_scbs)
}

# ==== scmet ====

summarizeResScmet <- function(NV) {
  
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
  
  smr_scmet <- data.frame(method = "scMET",
                          NV = NV,
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
  return(smr_scmet)
}

# ==== summarized plots ====
NV <- 5000
smr_vseq <- rbind(summarizeResVseq(NV))
smr_vseq_cr <- rbind(summarizeResVseqCr(NV))
smr_scbs <- rbind(summarizeResScbs(NV))
smr_scmet <- rbind(summarizeResScmet(NV))
smr <- rbind(smr_vseq[-3], smr_vseq_cr[-3], smr_scbs[-3], smr_scmet[-3])

colors <- RColorBrewer::brewer.pal(n = 4, name = "RdYlBu")

# site-level fdr and power
smr %>%  
  ggplot(aes(fdr_site, power_site, color = method)) + 
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dotted") +
  geom_path() +
  scale_color_manual(values = c("vmrseq" = colors[1], 
                                "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3],
                                "scMET" = colors[4])) +
  xlab("Site-level FDR") + ylab("Site-level power") +
  ggtitle("Modified real chromosome (200 cells, 4 subpops)") +
  xlim(0,1) + ylim(0,1)
ggsave(paste0(
  "plots/sim_studies/benchmark_real_chr/comparison/pseudoChr_fdr&power_siteLevel_", 
  N , "cells_", NP, "subpops_", 
  NV, "VMRs_seed", seed, ".png"
), width = 6, height = 5)

# region-level fdr and power
smr %>%  
  ggplot(aes(fdr_region, power_region, color = method)) + 
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dotted") +
  geom_path() +
  scale_color_manual(values = c("vmrseq" = colors[1], 
                                "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3],
                                "scMET" = colors[4])) +
  xlab("Region-level FDR") + ylab("Region-level power") +
  ggtitle("Modified real chromosome (200 cells, 4 subpops)") +
  xlim(0,1) + ylim(0,1)
ggsave(paste0(
  "plots/sim_studies/benchmark_real_chr/comparison/pseudoChr_fdr&power_regionLevel_", 
  N , "cells_", NP, "subpops_", 
  NV, "VMRs_seed", seed, ".png"
), width = 6, height = 5)
