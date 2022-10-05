.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15/")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
devtools::load_all("../vmrseq-package/vmrseq/")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(GenomicRanges))

NV <- 2000
seed <- 2022
chromosome <- "chr1"
subtype <- "IT-L23_Cux1"
min_olap <- 3 # minimum number of sites that overlaps with true VMR to be counted as true positive

# ==== vmrseq ====

.vmrseqEval <- function(res_vseq) {
  res_gr <- res_vseq$gr
  power_site <- with(res_gr, sum(is_vml&!is.na(vmr_index))/sum(is_vml))
  fdr_site <- with(res_gr, sum(!is_vml&!is.na(vmr_index))/sum(!is.na(vmr_index)))
  total_n_site <- with(res_gr, sum(!is.na(vmr_index)))

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
  power_region <- sum(true$n_dectected >= min_olap) / nrow(true)
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
  fdr_region <- 1 - sum(detected$n_true >= min_olap) / nrow(detected)
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
  power_region <- sum(true$n_dectected >= min_olap) / nrow(true)
  
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
  fdr_region <- 1 - sum(detected$n_true >= min_olap) / nrow(detected)
  return(list(total_n_site = total_n_site,
              power_site = power_site, fdr_site = fdr_site, 
              power_region = power_region , fdr_region  = fdr_region , 
              true = true, detected = detected))
}

summarizeResVseq <- function(NV, sparseLevel) {
  
  smr_vseq <- data.frame(method = "vmrseq",
                         NV = NV,
                         threshold = 0.05,
                         # threshold = c(seq(0.002, 0.005, 0.001), seq(0.01, 0.1, 0.01), 
                         #               0.12, 0.15, seq(0.2,0.4,0.1)),
                         total_n_site = NA,
                         power_site = NA,
                         fdr_site = NA,
                         power_region = NA,
                         fdr_region = NA)
  
  for (i in 1:nrow(smr_vseq)) {
    res_vseq <- readRDS(
      paste0(
        "data/interim/sim_studies/benchmark_real_chr/vmrseq/output/pseudoChr_",
        subtype, "_", chromosome, "_", 
        N, "cells_", NP, "subpops_", 
        NV, "VMRs_sparseLevel", sparseLevel,            
        "_alpha", smr_vseq$threshold[i], "_seed", seed, "_vmrseqOutput.rds"
      )
    )
    if (!is.null(res_vseq)) {
      smr <- .vmrseqEval(res_vseq)
      # View(smr$true); View(smr$detected)
      smr_vseq$total_n_site[i] <- smr$total_n_site
      smr_vseq$power_site[i] <- smr$power_site
      smr_vseq$fdr_site[i] <- smr$fdr_site
      smr_vseq$power_region[i] <- smr$power_region
      smr_vseq$fdr_region[i] <- smr$fdr_region
    }
    # cat(i, " ")
  }
  
  return(smr_vseq)
}

summarizeResVseqCr <- function(NV, sparseLevel) {
  smr_vseq_cr <- data.frame(method = "vmrseq_CR",
                            NV = NV,
                            threshold = 0.05,
                            # threshold = c(seq(0.002, 0.005, 0.001), seq(0.01, 0.1, 0.01), 
                            #               0.12, 0.15, seq(0.2,0.4,0.1)),
                            total_n_site = NA,
                            power_site = NA,
                            fdr_site = NA,
                            power_region = NA,
                            fdr_region = NA)
  
  for (i in 1:nrow(smr_vseq_cr)) {
    res_vseq <- readRDS(
      paste0(
        "data/interim/sim_studies/benchmark_real_chr/vmrseq/output/pseudoChr_",
        subtype, "_", chromosome, "_", 
        N, "cells_", NP, "subpops_", 
        NV, "VMRs_sparseLevel", sparseLevel, 
        "_alpha", smr_vseq$threshold[i], "_seed", seed, "_vmrseqOutput.rds"
      )
    )
    if (!is.null(res_vseq)) {
      smr <- .vmrseqCrEval(res_vseq)
      smr_vseq_cr$total_n_site[i] <- smr$total_n_site
      smr_vseq_cr$power_site[i] <- smr$power_site
      smr_vseq_cr$fdr_site[i] <- smr$fdr_site
      smr_vseq_cr$power_region[i] <- smr$power_region
      smr_vseq_cr$fdr_region[i] <- smr$fdr_region
    }
    # cat(i, " ")
  }
  return(smr_vseq_cr)
}

# === smallwood ====

summarizeResSmwd <- function(NV, sparseLevel) {
  gr <- loadHDF5SummarizedExperiment(paste0(
    "data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
    subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
    NV, "VMRs_sparseLevel", sparseLevel,            
    "_seed", seed)) %>% granges()
  res_smwd.gr <- readRDS(paste0(
    "data/interim/sim_studies/benchmark_real_chr/smallwood/output/pseudoChr_",
    subtype, "_", chromosome, "_",
    N, "cells_", NP, "subpops_",
    NV, "VMRs_sparseLevel", sparseLevel,            
    "_seed", seed, "_varLowerBound.rds")) 
  
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
    smr_smwd$power_region[i] <- sum(true_olap >= min_olap)/length(true_olap)
    
    detected_olap <- gr %>% as.data.frame %>% group_by(vmr_index) %>%  
      summarise(n_olap = sum(!is.na(vmr_name))) %>% dplyr::select(n_olap) %>% unlist
    smr_smwd$fdr_region[i] <- sum(detected_olap < min_olap)/length(detected_olap)
    
    # cat(i, " ")
  }
  return(smr_smwd)
}

# === scbs ====

summarizeResScbs <- function(NV, sparseLevel) {
  gr <- loadHDF5SummarizedExperiment(
    paste0(
      "data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
      subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
      NV, "VMRs_sparseLevel", sparseLevel,            
      "_seed", seed
    )
  ) %>% granges()
  
  smr_scbs <- data.frame(method = "scbs",
                         NV = NV,
                         threshold = c(0.001, 0.005, seq(0.01, 0.025, 0.005), seq(.03, .09, .01), seq(0.1,0.3,0.1)),
                         total_n_site = NA,
                         power_site = NA,
                         fdr_site = NA,
                         power_region = NA,
                         fdr_region = NA)
  
  for (i in 1:nrow(smr_scbs)) {
    res_scbs <- fread(
      paste0("data/interim/sim_studies/benchmark_real_chr/scbs/output/pseudoChr_IT-L23_Cux1_chr1_", 
             N, "cells_", NP, "subpops_", 
             NV, "VMRs_sparseLevel", sparseLevel,            
             "_seed", seed, "_", 
             format(smr_scbs$threshold[i], scientific = FALSE), 
             "vt.bed")
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
    smr_scbs$power_region[i] <- sum(true_olap >= min_olap)/length(true_olap)
    
    detected_olap <- gr %>% as.data.frame %>% group_by(vmr_index) %>%  
      summarise(n_olap = sum(!is.na(vmr_name))) %>% dplyr::select(n_olap) %>% unlist
    smr_scbs$fdr_region[i] <- sum(detected_olap < min_olap)/length(detected_olap)
    
    
    # cat(i, " ")
  }
  return(smr_scbs)
}

# ==== scmet ====

summarizeResScmet <- function(NV, sparseLevel) {
  
  bp_size <- 20000
  gr <- loadHDF5SummarizedExperiment(
    paste0(
      "data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
      subtype, "_", chromosome, "_", N, "cells_", NP, "subpops_", 
      NV, "VMRs_sparseLevel", sparseLevel,            
      "_seed", seed
    )
  ) %>% granges()
  wds.gr <- readRDS(paste0(
    "data/interim/sim_studies/benchmark_real_chr/scmet/input/features_", 
    subtype, "_", chromosome, "_", 
    bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
    NV, "VMRs_sparseLevel", sparseLevel,            
    "_seed", seed, ".rds"
  ))
  
  smr_scmet <- data.frame(method = "scmet",
                          NV = NV,
                          threshold = c(0.02, 0.05, seq(0.1,0.9,0.1), 0.99),
                          total_n_site = NA,
                          power_site = NA,
                          fdr_site = NA,
                          power_region = NA,
                          fdr_region = NA)
  for (i in 1:nrow(smr_scmet)) {
    
    efdr <- smr_scmet$threshold[i]
    res_scmet <- fread(paste0(
      "data/interim/sim_studies/benchmark_real_chr/scmet/output/pseudoChr_",
      subtype, "_", chromosome, "_", 
      bp_size/1000, "kbWindow_", N, "cells_", NP, "subpops_", 
      NV, "VMRs_sparseLevel", sparseLevel,            
      "_seed", seed, "_", efdr, "efdr_scmetOutput.csv"
    )) %>% mutate(window = as.integer(gsub(".*_(.*)$", "\\1", feature_name))) %>% filter(is_variable)
    
    hvf.gr <- wds.gr[res_scmet$window]
    
    hits <- findOverlaps(gr, hvf.gr)
    gr$vmr_index <- NA # to store detected VMRs by scbs
    gr$vmr_index[queryHits(hits)] <- subjectHits(hits) 
    smr_scmet$power_site[i] <- with(gr, sum(is_vml&!is.na(vmr_index)) / sum(is_vml))
    smr_scmet$fdr_site[i] <- with(gr, sum(!is_vml&!is.na(vmr_index)) / sum(!is.na(vmr_index)))
    smr_scmet$total_n_site[i] <- with(gr, sum(!is.na(vmr_index)))
    
    true_olap <- gr %>% as.data.frame %>% group_by(vmr_name) %>%  
      summarise(n_olap = sum(!is.na(vmr_index))) %>% select(n_olap) %>% unlist
    smr_scmet$power_region[i] <- sum(true_olap >= min_olap)/length(true_olap)
    
    detected_olap <- gr %>% as.data.frame %>% group_by(vmr_index) %>%  
      summarise(n_olap = sum(!is.na(vmr_name))) %>% select(n_olap) %>% unlist
    smr_scmet$fdr_region[i] <- sum(detected_olap < min_olap)/length(detected_olap)
  }
  return(smr_scmet)
}

# === (Archived) Smooth vs. no smooth in vmrseq ====
# NV <- 2000
# smr1 <- summarizeResVseq(NV, penalty = 0, smooth = T) %>% mutate(smoothed = "TRUE")
# smr1_cr <- summarizeResVseqCr(NV, smooth = T)  %>% mutate(smoothed = "TRUE")
# smr0 <- summarizeResVseq(NV, penalty = 0, smooth = F) %>% mutate(smoothed = "FALSE")
# smr0_cr <- summarizeResVseqCr(NV, smooth = F)  %>% mutate(smoothed = "FALSE")
# smr <- rbind(smr1, smr1_cr, smr0, smr0_cr)
# 
# # site-level fdr and power
# smr %>%  
#   ggplot(aes(fdr_site, power_site, linetype = method, color = smoothed)) + 
#   geom_vline(xintercept = 0.05, linetype = "dotted") +
#   geom_path() + geom_point() +
#   xlab("Site-level FDR") + ylab("Site-level power") +
#   ggtitle(paste0("Modified real chromosome (", N, " cells, ", NP, " subpops)")) +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
#   theme_classic()
# ggsave(paste0(
#   "plots/sim_studies/benchmark_real_chr/comparison/ifSmooth_fdr&power_siteLevel_", 
#   N , "cells_", NP, "subpops_", 
#   NV, "VMRs_seed", seed, ".png"
# ), width = 8, height = 6)
# 
# # region-level fdr and power
# smr %>%  
#   ggplot(aes(fdr_region, power_region, linetype = method, color = smoothed)) + 
#   geom_vline(xintercept = 0.05, linetype = "dotted") +
#   geom_path() + geom_point() +
#   xlab("Region-level FDR") + ylab("Region-level power") +
#   ggtitle(paste0("Modified real chromosome (", N, " cells, ", NP, " subpops)")) +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
#   theme_classic()
# ggsave(paste0(
#   "plots/sim_studies/benchmark_real_chr/comparison/ifSmooth_fdr&power_regionLevel_", 
#   N , "cells_", NP, "subpops_", 
#   NV, "VMRs_seed", seed, ".png"
# ), width = 8, height = 6)

# ==== Compare methods ====
N <- 500
for (NP in c(2,3,4,5,8,12,20)){
  for (sparseLevel in 1:3) {
    smr_vseq <- summarizeResVseq(NV, sparseLevel)
    smr_vseq_cr <- summarizeResVseqCr(NV, sparseLevel) 
    smr_scbs <- summarizeResScbs(NV, sparseLevel)
    smr_smwd <- summarizeResSmwd(NV, sparseLevel)
    smr_scmet <- summarizeResScmet(NV, sparseLevel)
    smr <- rbind(smr_vseq, smr_vseq_cr, smr_scbs, smr_smwd, smr_scmet)
    
    fwrite(smr, paste0(
      "data/interim/sim_studies/benchmark_real_chr/summary/pseudoChr_fdr&power_siteLevel_",
      N , "cells_", NP, "subpops_",
      NV, "VMRs_sparseLevel", sparseLevel,
      "_seed", seed, ".csv.gz"
    ))
    
    colors <- RColorBrewer::brewer.pal(n = 6, name = "RdYlBu")[-4]
    
    # site-level fdr and power
    smr %>%  
      ggplot(aes(fdr_site, power_site, color = method)) + 
      geom_vline(xintercept = 0.05, linetype = "dotted") +
      geom_path() + geom_point() +
      scale_color_manual(values = c("vmrseq" = colors[1], 
                                    "vmrseq_CR" = colors[2], 
                                    "smallwood" = colors[3],
                                    "scbs" = colors[4],
                                    "scmet" = colors[5])) +
      xlab("Site-level FDR") + ylab("Site-level power") +
      ggtitle(paste0("Modified real chromosome (",N," cells, ",NP," subpops, sparse level ",sparseLevel,")")) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
      theme_classic()
    ggsave(paste0(
      "plots/sim_studies/benchmark_real_chr/comparison/pseudoChr_fdr&power_siteLevel_",
      N , "cells_", NP, "subpops_",
      NV, "VMRs_sparseLevel", sparseLevel,
      "_seed", seed, ".png"
    ), width = 8, height = 7.5)
    
    # region-level fdr and power
    smr %>%  
      ggplot(aes(fdr_region, power_region, color = method)) + 
      geom_vline(xintercept = 0.05, linetype = "dotted") +
      geom_path() + geom_point() +
      scale_color_manual(values = c("vmrseq" = colors[1], 
                                    "vmrseq_CR" = colors[2], 
                                    "smallwood" = colors[3],
                                    "scbs" = colors[4],
                                    "scmet" = colors[5])) +
      xlab("Region-level FDR") + ylab("Region-level power") +
      ggtitle(paste0("Modified real chromosome (",N," cells, ",NP," subpops, sparse level ",sparseLevel,")")) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
      theme_classic()
    ggsave(paste0(
      "plots/sim_studies/benchmark_real_chr/comparison/pseudoChr_fdr&power_regionLevel_",
      N , "cells_", NP, "subpops_",
      NV, "VMRs_sparseLevel", sparseLevel,
      "_seed", seed, ".png"
    ), width = 8, height = 7.5)
    print(paste0("NP = ", NP, ", sparseLevel = ", sparseLevel))
  }
}
