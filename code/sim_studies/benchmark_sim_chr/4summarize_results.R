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
subtype <- "IT-L23_Cux1"

# ==== vmrseq ====
.vmrseqEval <- function(res_vseq) {
  res_gr <- res_vseq$gr
  # index <- which(res_gr$loglik_diff < penalty) # indices of sites failed to pass penalty threshold
  # res_gr$vmr_index[index] <- NA
  # res_gr$vmr_num_cpg[index] <- NA
  power_site <- with(res_gr, sum(is_vml&!is.na(vmr_index))/sum(is_vml))
  fdr_site <- with(res_gr, sum(!is_vml&!is.na(vmr_index))/sum(!is.na(vmr_index)))
  total_n_site <- with(res_gr, sum(!is.na(vmr_index)))
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
              med_total = median(total)) %>%
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
              mean_MF = mean(meth/total))
  fdr_region <- 1-sum(detected$n_true > 0) / nrow(detected)
  # cat("Region-level FDR =", fdr, "\n")
  return(list(total_n_site = total_n_site,
              power_site = power_site, fdr_site = fdr_site, 
              power_region = power_region , fdr_region  = fdr_region , 
              true = true, detected = detected))
}

.vmrseqCrEval <- function(res_vseq) {
  res_gr <- res_vseq$gr
  power_site <- with(res_gr, sum(is_vml&!is.na(cr_index))/sum(is_vml))
  fdr_site <- with(res_gr, sum(!is_vml&!is.na(cr_index))/sum(!is.na(cr_index)))
  total_n_site <- with(res_gr, sum(!is.na(cr_index)))
  
  true <- res_gr %>% data.frame %>%
    filter(is_vml) %>%
    group_by(vmr_name) %>%
    summarise(n_cpg = n(),
              n_dectected = sum(!is.na(cr_index)),
              n_incr = sum(!is.na(cr_index)),
              cr_index = max(cr_index, na.rm = T),
              pi1 = max(pi1, na.rm = T),
              med_total = median(total)) %>%
    filter(n_cpg >= 5)
  power_region <- sum(true$n_dectected > 0) / nrow(true)
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
              mean_MF = mean(meth/total))
  fdr_region <- 1-sum(detected$n_true > 0) / nrow(detected)
  return(list(total_n_site = total_n_site,
              power_site = power_site, fdr_site = fdr_site, 
              power_region = power_region , fdr_region  = fdr_region , 
              true = true, detected = detected))
}

summarizeResVseq <- function(NV) {
  
  smr_vseq <- data.frame(method = "vmrseq",
                         NV = NV,
                         qVar = c(0.005, 
                                  seq(0.01, 0.1, 0.01), 
                                  0.12, 0.15, seq(0.2,0.4,0.1)),
                         total_n_site = NA,
                         power_site = NA,
                         fdr_site = NA,
                         power_region = NA,
                         fdr_region = NA)
  
  for (i in 1:nrow(smr_vseq)) {
    res_vseq <- readRDS(
      paste0(
        "data/interim/sim_studies/benchmark_sim_chr/vmrseq/output/simChr_",
        subtype, "_", chromosome, "_", 
        N, "cells_", NP, "subpops_", 
        NV, "VMRs_qVar", smr_vseq$qVar[i], "_seed", seed, "_vmrseqOutput.rds"
      )
    )
    smr <- .vmrseqEval(res_vseq)
    # View(smr$true); View(smr$detected)
    smr_vseq$total_n_site[i] <- smr$total_n_site
    smr_vseq$power_site[i] <- smr$power_site
    smr_vseq$fdr_site[i] <- smr$fdr_site
    smr_vseq$power_region[i] <- smr$power_region
    smr_vseq$fdr_region[i] <- smr$fdr_region
    cat(i, " ")
  }
  
  return(smr_vseq)
}

summarizeResVseqCr <- function(NV) {
  smr_vseq_cr <- data.frame(method = "vmrseq_CR",
                            NV = NV,
                            qVar = c(0.005, 
                                     seq(0.01, 0.1, 0.01), 
                                     0.12, 0.15, seq(0.2,0.4,0.1)),
                            total_n_site = NA,
                            power_site = NA,
                            fdr_site = NA,
                            power_region = NA,
                            fdr_region = NA)
  
  for (i in 1:nrow(smr_vseq_cr)) {
    res_vseq <- readRDS(
      paste0(
        "data/interim/sim_studies/benchmark_sim_chr/vmrseq/output/simChr_",
        subtype, "_", chromosome, "_", 
        N, "cells_", NP, "subpops_", 
        NV, "VMRs_qVar", smr_vseq_cr$qVar[i], "_seed", seed, "_vmrseqOutput.rds"
      )
    )
    smr <- .vmrseqCrEval(res_vseq)
    smr_vseq_cr$total_n_site[i] <- smr$total_n_site
    smr_vseq_cr$power_site[i] <- smr$power_site
    smr_vseq_cr$fdr_site[i] <- smr$fdr_site
    smr_vseq_cr$power_region[i] <- smr$power_region
    smr_vseq_cr$fdr_region[i] <- smr$fdr_region
    cat(i, " ")
  }
  return(smr_vseq_cr)
}

# === smallwood ====
summarizeResSmwd <- function(NV) {
  gr <- loadHDF5SummarizedExperiment(paste0(
    "data/interim/sim_studies/benchmark_sim_chr/simulated/simChr_",
    subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
    NV, "VMRs_seed", seed)) %>% granges()
  res_smwd.gr <- readRDS(paste0(
    "data/interim/sim_studies/benchmark_sim_chr/smallwood/output/simChr_",
    subtype, "_", chromosome, "_",
    N, "cells_", NP, "subpops_",
    NV, "VMRs_seed", seed, "_varLowerBound.rds")) 
  
  smr_smwd <- data.frame(method = "smallwood",
                         NV = NV,
                         threshold = c(1e-4, 5e-4, 1e-3, 5e-3, seq(0.01, 0.025, 0.005), seq(.03, .09, .01), seq(0.1,0.3,0.1)),
                         total_n_site = NA,
                         power_site = NA,
                         fdr_site = NA,
                         power_region = NA,
                         fdr_region = NA)
  
  for (i in 1:nrow(smr_smwd)) {
    top_n <- round(length(res_smwd.gr)*smr_smwd$threshold[i])
    top_idx <- order(res_smwd.gr$var_lb, decreasing = TRUE)[1:top_n] %>% sort
    hvf.gr <- res_smwd.gr[top_idx]
    hits <- findOverlaps(gr, hvf.gr)
    gr$vmr_index <- NA # to store detected VMRs by scbs
    gr$vmr_index[queryHits(hits)] <- subjectHits(hits) 
    smr_smwd$power_site[i] <- with(gr, sum(is_vml&!is.na(vmr_index)) / sum(is_vml))
    smr_smwd$fdr_site[i] <- with(gr, sum(!is_vml&!is.na(vmr_index)) / sum(!is.na(vmr_index)))
    smr_smwd$total_n_site[i] <- with(gr, sum(!is.na(vmr_index)))
    
    true_olap <- gr %>% as.data.frame %>% group_by(vmr_name) %>%  
      summarise(n_olap = sum(!is.na(vmr_index))) %>% dplyr::select(n_olap) %>% unlist
    smr_smwd$power_region[i] <- sum(true_olap>0)/length(true_olap)
    
    detected_olap <- gr %>% as.data.frame %>% group_by(vmr_index) %>%  
      summarise(n_olap = sum(!is.na(vmr_name))) %>% dplyr::select(n_olap) %>% unlist
    smr_smwd$fdr_region[i] <- sum(detected_olap==0)/length(detected_olap)
    
    cat(i, " ")
  }
  return(smr_smwd)
}

# === scbs ====
summarizeResScbs <- function(NV) {
  gr <- loadHDF5SummarizedExperiment(
    paste0(
      "data/interim/sim_studies/benchmark_sim_chr/simulated/simChr_",
      subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
      NV, "VMRs_seed", seed
    )
  ) %>% granges()
  
  smr_scbs <- data.frame(method = "scbs",
                         NV = NV,
                         vt = c(0.001, 0.005, seq(0.01, 0.025, 0.005), seq(.03, .09, .01), seq(0.1,0.3,0.1)),
                         total_n_site = NA,
                         power_site = NA,
                         fdr_site = NA,
                         power_region = NA,
                         fdr_region = NA)
  
  for (i in 1:nrow(smr_scbs)) {
    res_scbs <- fread(
      paste0("data/interim/sim_studies/benchmark_sim_chr/scbs/output/simChr_IT-L23_Cux1_chr1_", 
             N, "cells_", NP, "subpops_", 
             NV, "VMRs_seed", seed, "_", format(smr_scbs$vt[i], scientific = FALSE), "vt.bed")
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
    smr_scbs$total_n_site[i] <- with(gr, sum(!is.na(vmr_index)))
    
    true_olap <- gr %>% as.data.frame %>% group_by(vmr_name) %>%  
      summarise(n_olap = sum(!is.na(vmr_index))) %>% dplyr::select(n_olap) %>% unlist
    smr_scbs$power_region[i] <- sum(true_olap>0)/length(true_olap)
    
    detected_olap <- gr %>% as.data.frame %>% group_by(vmr_index) %>%  
      summarise(n_olap = sum(!is.na(vmr_name))) %>% dplyr::select(n_olap) %>% unlist
    smr_scbs$fdr_region[i] <- sum(detected_olap==0)/length(detected_olap)
    
    
    cat(i, " ")
  }
  return(smr_scbs)
}

# ==== scmet ====
# summarizeResScmet <- function(NV) {
#   
#   bp_size <- 20000
#   gr <- loadHDF5SummarizedExperiment(
#     paste0(
#       "data/interim/sim_studies/benchmark_sim_chr/simulated/simChr_",
#       subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
#       NV, "VMRs_seed", seed
#     )
#   ) %>% granges()
#   wds.gr <- readRDS(paste0(
#     "data/interim/sim_studies/benchmark_sim_chr/scmet/input/features_", 
#     subtype, "_", chromosome, "_", 
#     bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
#     NV, "VMRs_seed", seed, ".rds"
#   ))
#   
#   smr_scmet <- data.frame(method = "scmet",
#                           NV = NV,
#                           efdr = c(0.02, 0.05, seq(0.1,0.9,0.1), 0.99),
#                           total_n_site = NA,
#                           power_site = NA,
#                           fdr_site = NA,
#                           power_region = NA,
#                           fdr_region = NA)
#   for (i in 1:nrow(smr_scmet)) {
#     
#     efdr <- smr_scmet$efdr[i]
#     res_scmet <- fread(paste0(
#       "data/interim/sim_studies/benchmark_sim_chr/scmet/output/simChr_",
#       subtype, "_", chromosome, "_", 
#       bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
#       NV, "VMRs_seed", seed, "_", efdr, "efdr_scmetOutput.csv"
#     )) %>% mutate(window = as.integer(gsub(".*_(.*)$", "\\1", feature_name))) %>% filter(is_variable)
#     
#     hvf.gr <- wds.gr[res_scmet$window]
#     
#     hits <- findOverlaps(gr, hvf.gr)
#     gr$vmr_index <- NA # to store detected VMRs by scbs
#     gr$vmr_index[queryHits(hits)] <- subjectHits(hits) 
#     smr_scmet$power_site[i] <- with(gr, sum(is_vml&!is.na(vmr_index)) / sum(is_vml))
#     smr_scmet$fdr_site[i] <- with(gr, sum(!is_vml&!is.na(vmr_index)) / sum(!is.na(vmr_index)))
#     smr_scmet$total_n_site[i] <- with(gr, sum(!is.na(vmr_index)))
#     
#     true_olap <- gr %>% as.data.frame %>% group_by(vmr_name) %>%  
#       summarise(n_olap = sum(!is.na(vmr_index))) %>% select(n_olap) %>% unlist
#     smr_scmet$power_region[i] <- sum(true_olap>0)/length(true_olap)
#     
#     detected_olap <- gr %>% as.data.frame %>% group_by(vmr_index) %>%  
#       summarise(n_olap = sum(!is.na(vmr_name))) %>% select(n_olap) %>% unlist
#     smr_scmet$fdr_region[i] <- sum(detected_olap==0)/length(detected_olap)
#   }
#   return(smr_scmet)
# }

# ==== Compare methods ====
NV <- 2000
smr_vseq <- summarizeResVseq(NV)
smr_vseq_cr <- summarizeResVseqCr(NV) 
smr_scbs <- summarizeResScbs(NV)
smr_smwd <- summarizeResSmwd(NV)
# smr_scmet <- summarizeResScmet(NV)
smr <- rbind(smr_vseq[-3], smr_vseq_cr[-3], smr_scbs[-3], smr_smwd[-3])
# smr <- rbind(smr_vseq[-3], smr_vseq_cr[-3], smr_scbs[-3])

colors <- RColorBrewer::brewer.pal(n = 6, name = "RdYlBu")[-4]

# site-level fdr and power
smr %>%  
  ggplot(aes(fdr_site, power_site, color = method)) + 
  geom_vline(xintercept = 0.05, linetype = "dotted") +
  geom_path() + geom_point() +
  scale_color_manual(values = c("vmrseq" = colors[1], 
                                "vmrseq_CR" = colors[2], 
                                "smallwood" = colors[3],
                                "scbs" = colors[4])) +
  xlab("Site-level FDR") + ylab("Site-level power") +
  ggtitle(paste0("Modified real chromosome (", N, " cells, ", NP, " subpops)")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  theme_classic()
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/simChr_fdr&power_siteLevel_", 
  N , "cells_", NP, "subpops_", 
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 7.5)

# region-level fdr and power
smr %>%  
  ggplot(aes(fdr_region, power_region, color = method)) + 
  geom_vline(xintercept = 0.05, linetype = "dotted") +
  geom_path() + geom_point() +
  scale_color_manual(values = c("vmrseq" = colors[1], 
                                "vmrseq_CR" = colors[2], 
                                "smallwood" = colors[3],
                                "scbs" = colors[4])) +
  xlab("Region-level FDR") + ylab("Region-level power") +
  ggtitle(paste0("Modified real chromosome (", N, " cells, ", NP, " subpops)")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  theme_classic()
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/simChr_fdr&power_regionLevel_", 
  N , "cells_", NP, "subpops_", 
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 7.5)

