source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

N <- 200

NV <- 2000
seed <- 2022
chromosome <- "chr1"
subtype <- "IT-L23_Cux1"
NPs <- c(2,3,4,5,8,12,20)

# colors <- RColorBrewer::brewer.pal(n = 6, name = "RdYlBu")[-4]
colors <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]

# ==== plot fdr&power at default threshold ====

getResultsDefaultSetting <- function(N) {
  
  getOneRow <- function(N, NP, sparseLevel) {
    smr0 <- fread(paste0(
      "data/interim/sim_studies/benchmark_sim_chr/summary/simChr_fdr&power_siteLevel_",
      N , "cells_", NP, "subpops_",
      NV, "VMRs_sparseLevel", sparseLevel,
      "_seed", seed, ".csv.gz"
    ))
    smr0_dft <- smr0 %>% filter(
      (method == "vmrseq" & threshold == 0.05) | 
        (method == "vmrseq_CR" & threshold == 0.05) | 
        (method == "scbs" & threshold == 0.02) | 
        (method == "smallwood" & threshold == 0.02)
    )
    return(data.frame("N" = N, "NP" = NP, "sparseLevel" = sparseLevel, smr0_dft))
  }
  
  smr_dft <- rbind(
    do.call(rbind, lapply(NPs, getOneRow, N = N, sparseLevel = 1)),
    do.call(rbind, lapply(NPs, getOneRow, N = N, sparseLevel = 2)),
    do.call(rbind, lapply(NPs, getOneRow, N = N, sparseLevel = 3))
  )  %>% 
    mutate(sparseLevel = factor(sparseLevel),
           method = factor(method, levels = c("vmrseq", "vmrseq_CR", "scbs", "smallwood")),
           F1_like_site = 2/(1/power_site + 1/(1-fdr_site)),
           F1_like_region = 2/(1/power_region + 1/(1-fdr_region)))
  return(smr_dft)
}

smr_dft <- getResultsDefaultSetting(N = N)

# site-level F1 score
smr_dft %>% ggplot(aes(color = method)) +
  geom_point(aes(NP, F1_like_site)) +
  geom_line(aes(NP, F1_like_site)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + scale_x_continuous(trans = "log2", breaks = NPs) +
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + ylab("Site-level F1-like score") + xlab("Number of subpopulations") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/point_siteLevel_f1_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)
smr_dft %>% ggplot(aes(color = method)) +
  geom_boxplot(aes(method, F1_like_site)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + 
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Site-level F1-like score") + xlab("Sparsity level") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/boxplot_siteLevel_f1_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)

# site-level FDR
smr_dft %>% ggplot(aes(color = method)) +
  geom_point(aes(NP, fdr_site)) +
  geom_line(aes(NP, fdr_site)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + scale_x_continuous(trans = "log2", breaks = NPs) +
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + ylab("Site-level FDR") + xlab("Number of subpopulations") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/point_siteLevel_fdr_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)
smr_dft %>% ggplot(aes(color = method)) +
  geom_boxplot(aes(method, fdr_site)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) +
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Site-level FDR") + xlab("Sparsity level") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/boxplot_siteLevel_fdr_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)

# site-level power
smr_dft %>% ggplot(aes(color = method)) +
  geom_point(aes(NP, power_site)) +
  geom_line(aes(NP, power_site)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + scale_x_continuous(trans = "log2", breaks = NPs) +
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + ylab("Site-level power") + xlab("Number of subpopulations") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/point_siteLevel_power_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)
smr_dft %>% ggplot(aes(color = method)) +
  geom_boxplot(aes(method, power_site)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + 
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Site-level power") + xlab("Number of subpopulations") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/boxplot_siteLevel_power_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)

# region-level F1 score
smr_dft %>% ggplot(aes(color = method)) +
  geom_point(aes(NP, F1_like_region)) +
  geom_line(aes(NP, F1_like_region)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + scale_x_continuous(trans = "log2", breaks = NPs) +
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + ylab("Region-level F1-like score") + xlab("Number of subpopulations") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/point_regionLevel_f1_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)
smr_dft %>% ggplot(aes(color = method)) +
  geom_boxplot(aes(method, F1_like_region)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + 
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Region-level F1-like score") + xlab("Sparsity level") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/boxplot_regionLevel_f1_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)

# region-level FDR
smr_dft %>% ggplot(aes(color = method)) +
  geom_point(aes(NP, fdr_region)) +
  geom_line(aes(NP, fdr_region)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + scale_x_continuous(trans = "log2", breaks = NPs) +
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + ylab("Region-level FDR") + xlab("Number of subpopulations") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/point_regionLevel_fdr_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)
smr_dft %>% ggplot(aes(color = method)) +
  geom_boxplot(aes(method, fdr_region)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + 
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Region-level FDR") + xlab("Sparsity level") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/boxplot_regionLevel_fdr_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)

# region-level power
smr_dft %>% ggplot(aes(color = method)) +
  geom_point(aes(NP, power_region)) +
  geom_line(aes(NP, power_region)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + scale_x_continuous(trans = "log2", breaks = NPs) +
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + ylab("Region-level power") + xlab("Number of subpopulations") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/point_regionLevel_power_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)
smr_dft %>% ggplot(aes(color = method)) +
  geom_boxplot(aes(method, power_region)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + 
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Region-level power") + xlab("Sparsity level") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/boxplot_regionLevel_power_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)


# ==== plot partial AUC (fdr <= 0.2) ====

computeRRA <- function(power, fdr, xlim = c(0, 0.2)) { # The Ratio of Relevant Areas (RRA) indicator
  
  idx <- which(fdr >= xlim[1] & fdr <= xlim[2])
  pfdr <- c(xlim[1], fdr[idx], xlim[2])
  ppower <- c(power[idx], power[idx[length(idx)]])
  pauc <- sum(ppower * diff(pfdr))
  rra <- pauc / (xlim[2] - xlim[1])
  
  return(rra)
} 

getRRA <- function(N) {
  
  computeRRAbyExpr <- function(N, NP, sparseLevel) {
    
    smr0 <- fread(paste0(
      "data/interim/sim_studies/benchmark_sim_chr/summary/simChr_fdr&power_siteLevel_",
      N , "cells_", NP, "subpops_",
      NV, "VMRs_sparseLevel", sparseLevel,
      "_seed", seed, ".csv.gz"
    ))
    
    smr0_rra <- smr0 %>% na.omit() %>% 
      group_by(method) %>%
      summarise(rra_site = computeRRA(power_site, fdr_site),
                rra_region = computeRRA(power_region, fdr_region))
    return(data.frame("N" = N, "NP" = NP, "sparseLevel" = sparseLevel, smr0_rra))
  }
  
  smr_rra <- rbind(
    do.call(rbind, lapply(NPs, computeRRAbyExpr, N = N, sparseLevel = 1)),
    do.call(rbind, lapply(NPs, computeRRAbyExpr, N = N, sparseLevel = 2)),
    do.call(rbind, lapply(NPs, computeRRAbyExpr, N = N, sparseLevel = 3))
  )  %>% 
    mutate(sparseLevel = factor(sparseLevel),
           method = factor(method, levels = c("vmrseq", "vmrseq_CR", "scbs", "smallwood")))
  return(smr_rra)
}

smr_rra <- getRRA(N = N)

# site-level RRA indicator
smr_rra %>% ggplot(aes(color = method)) +
  geom_point(aes(NP, rra_site)) +
  geom_line(aes(NP, rra_site)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + scale_x_continuous(trans = "log2", breaks = NPs) +
  ggtitle(paste0("Simulated chromosome (",N," cells)")) + 
  theme_classic() + ylab("Site-level RRA") + xlab("Number of subpopulations") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/point_siteLevel_RRA_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)
smr_rra %>% ggplot(aes(color = method)) +
  geom_boxplot(aes(method, rra_site)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Site-level RRA") + xlab("Sparsity level") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/boxplot_siteLevel_RRA_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)

# region-level RRA indicator
smr_rra %>% ggplot(aes(color = method)) +
  geom_point(aes(NP, rra_region)) +
  geom_line(aes(NP, rra_region)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + scale_x_continuous(trans = "log2", breaks = NPs) +
  ggtitle(paste0("Simulated chromosome (",N," cells)")) + 
  theme_classic() + ylab("Region-level RRA") + xlab("Number of subpopulations") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/point_regionLevel_RRA_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)
smr_rra %>% ggplot(aes(color = method)) +
  geom_boxplot(aes(method, rra_region)) +
  facet_wrap(~ sparseLevel) +
  scale_color_manual(values = c("vmrseq" = colors[1], "vmrseq_CR" = colors[2], 
                                "scbs" = colors[3], "smallwood" = colors[4])) + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  ggtitle(paste0("Simulated chromosome (",N," cells)")) +
  theme_classic() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ylab("Region-level RRA") + xlab("Sparsity level") 
ggsave(paste0(
  "plots/sim_studies/benchmark_sim_chr/comparison/boxplot_regionLevel_RRA_", N , "cells_",
  NV, "VMRs_seed", seed, ".png"
), width = 8, height = 4)





