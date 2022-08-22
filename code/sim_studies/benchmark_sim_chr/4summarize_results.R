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
  
  smr_true <- res_gr %>% data.frame %>%
    filter(is_vml) %>%
    group_by(vmr_name) %>%
    summarise(n_cpg = n(),
              n_dectected = sum(!is.na(vmr_index)),
              n_incr = sum(!is.na(cr_index)),
              cr_index = max(cr_index, na.rm = T),
              pi1 = max(pi1, na.rm = T),
              med_total = median(total),
              loglik_diff = max(loglik_diff, na.rm = T))
  power_region <- sum(smr_true$n_dectected > 0) / nrow(smr_true)
  # cat("Region-level power =", power, "; ")
  smr_detected <- res_gr %>% data.frame %>%
    filter(!is.na(vmr_index)) %>%
    group_by(vmr_index) %>%
    summarise(cr_index = max(cr_index),
              n_cpg = n(),
              n_true = sum(!is.na(vmr_name)),
              pi1 = max(pi1, na.rm = T),
              optim_pi = max(pi),
              mean_MF = mean(meth/total),
              loglik_diff = max(loglik_diff))
  fdr_region <- 1-sum(smr_detected$n_true > 0) / nrow(smr_detected)
  # cat("Region-level FDR =", fdr, "\n")
  return(list(power_site = power_site, fdr_site = fdr_site, 
              power_region = power_region , fdr_region  = fdr_region , 
              true = smr_true, detected = smr_detected))
}

res_vmrseq <- readRDS(
  paste0("data/interim/sim_studies/benchmark_sim_chr/vmrseq/output/simChr_IT-L23_Cux1_chr1_",
         N, "cells_", NP, "subpops_seed", seed, "_vmrseqOutput.rds")
)

smr_vmrseq <- data.frame(method = "vmrseq",
                         penalty = c(0:10, seq(15,60,5)),
                         power_site = NA,
                         fdr_site = NA,
                         power_region = NA,
                         fdr_region = NA)

for (i in 1:nrow(smr_vmrseq)) {
  pen <- smr_vmrseq$penalty[i]
  smr <- vmrseq_eval(res_vmrseq, penalty = pen)
  # if (pen==0) {View(smr$true); View(smr$detected)}
  smr_vmrseq$penalty[i] <- pen
  smr_vmrseq$power_site[i] <- smr$power_site
  smr_vmrseq$fdr_site[i] <- smr$fdr_site
  smr_vmrseq$power_region[i] <- smr$power_region
  smr_vmrseq$fdr_region[i] <- smr$fdr_region
  cat(i, " ")
}

# with(smr_vmrseq, plot(fdr_site, power_site, type = "l", xlim = c(0,1), ylim = c(0,1)))
# with(smr_vmrseq, plot(fdr_region, power_region, type = "l", xlim = c(0,1), ylim = c(0,1)))

# === scbs ====
gr <- loadHDF5SummarizedExperiment(
  paste0("data/interim/sim_studies/benchmark_sim_chr/simulated/simChr_IT-L23_Cux1_chr1_",
         N , "cells_", NP, "subpops_seed", seed)
) %>% granges()

smr_scbs <- data.frame(method = "scbs",
                       vt = c(seq(.001, .009, .001), seq(.01, .09, .01), seq(0.1,0.3,0.1)),
                       power_site = NA,
                       fdr_site = NA,
                       power_region = NA,
                       fdr_region = NA)

for (i in 1:nrow(smr_scbs)) {
  vt <- smr_scbs$vt[i]
  res_scbs <- fread(
    paste0("data/interim/sim_studies/benchmark_sim_chr/scbs/output/simChr_IT-L23_Cux1_chr1_", 
           N, "cells_", NP, "subpops_seed", seed, "_",vt,"vt.bed")
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

# res_scmet <- fread("data/interim/sim_studies/benchmark_sim_chr/scmet/output/")




# ==== plots ====
smr <- rbind(smr_vmrseq[-2], smr_scbs[-2])

# site-level fdr and power
smr %>%  
  ggplot(aes(fdr_site, power_site, color = method)) + 
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dashed") +
  geom_path() 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/simChr_fdr&power_siteLevel_", 
  N , "cells_", NP, "subpops_seed", seed, ".png"
))

# region-level fdr and power
smr %>%  
  ggplot(aes(fdr_region, power_region, color = method)) + 
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dashed") +
  geom_path() 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/simChr_fdr&power_regionLevel_", 
  N , "cells_", NP, "subpops_seed", seed, ".png"
))
