source("code/SETPATHS.R")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)

read_dir <- "data/interim/sim_studies/"
plot_dir <- "manuscript_related/manuscript_figures/simulation/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

# ---- fixed arguments ----
NV <- 2000
seed <- 2022
chromosome <- "chr1"
subtype <- "IT-L23_Cux1"
NPs <- c(2,3,4,5,8,12,20)

methodName <- function(method) switch (method,
                                       'vmrseq' = 'vmrseq',
                                       'vmrseq_CR' = 'vmrseq CRs',
                                       'scbs' = 'scbs',
                                       'smallwood' = 'Smallwood',
                                       'scmet' = 'scMET')
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = COLORS[1], "vmrseq CRs" = COLORS[2],
                 "scbs" = COLORS[3], "Smallwood" = COLORS[4], "scMET" = COLORS[5])
LEVELS <- c('vmrseq', 'vmrseq CRs', 'scbs', 'Smallwood', 'scMET')

# ---- Obtain legend ----

experiment <- 'real'; N = 200; NP = 2; sparseLevel = 1
smr0 <- fread(here(read_dir, 
                   paste0("benchmark_", experiment, "_chr"), 
                   "summary", 
                   paste0(ifelse(experiment == "real", "pseudo", "sim"),
                          "Chr_fdr&power_siteLevel_",
                          N , "cells_", NP, "subpops_",
                          NV, "VMRs_sparseLevel", sparseLevel,
                          "_seed", seed, ".csv.gz")))
p <- smr0 %>% 
  group_by(method) %>%
  summarise(min := min(fdr_site), med := median(fdr_site), max := max(fdr_site)) %>%
  ggplot(aes(method, med, color = method)) +
  geom_point() +
  geom_errorbar(aes(ymin = min, ymax = max)) +
  scale_color_manual(values = COLORVALUES) +
  theme_classic() + 
  theme(legend.position = 'bottom') +
  labs(color = 'Method')
legend <- cowplot::get_legend(p)
png(here(plot_dir, "legend_methods.png"), width = 2000, height = 200, res = 300)
grid::grid.newpage()
grid::grid.draw(legend)
dev.off()

# ---- F1-score, precision, and recall ----

plotDefaultSetting <- function(experiment) {
  
  getResultsDefaultSetting <- function(N) {
    getOneRow <- function(N, NP, sparseLevel) {
      smr0 <- fread(here(read_dir, 
                         paste0("benchmark_", experiment, "_chr"), 
                         "summary", 
                         paste0(ifelse(experiment == "real", "pseudo", "sim"),
                                "Chr_fdr&power_siteLevel_",
                                N , "cells_", NP, "subpops_",
                                NV, "VMRs_sparseLevel", sparseLevel,
                                "_seed", seed, ".csv.gz")))
      smr0_dft <- smr0 %>% filter(
        (method == "vmrseq" & threshold == 0.05) | 
          (method == "vmrseq_CR" & threshold == 0.05) | 
          (method == "scbs" & threshold == 0.02) | 
          (method == "smallwood" & threshold == 0.02) |
          (method == "scmet" & threshold == 0.1) 
      )
      return(data.frame("N" = N, "NP" = NP, "sparseLevel" = sparseLevel, smr0_dft))
    }
    
    smr_dft <- rbind(
      do.call(rbind, lapply(NPs, getOneRow, N = N, sparseLevel = 1)),
      do.call(rbind, lapply(NPs, getOneRow, N = N, sparseLevel = 2)),
      do.call(rbind, lapply(NPs, getOneRow, N = N, sparseLevel = 3))
    )  %>% 
      mutate(N = factor(paste0('N cells = ', N), levels = c('N cells = 200', 'N cells = 500', 'N cells = 1000', 'N cells = 2000')),
             sparseLevel = factor(sparseLevel, labels = c('high', 'medium', 'low')),
             method = factor(sapply(method, methodName), levels = LEVELS),
             F1_site = 2/(1/power_site + 1/(1-fdr_site)),
             F1_region = 2/(1/power_region + 1/(1-fdr_region)))
    return(smr_dft)
  }
  
  smr_dft0 <- rbind(
    getResultsDefaultSetting(N = 200),
    getResultsDefaultSetting(N = 500),
    getResultsDefaultSetting(N = 1000),
    getResultsDefaultSetting(N = 2000)
  ) 
  
  mySummarise <- function(col_name) {
    p <- smr_dft0 %>%
      group_by(N, sparseLevel, method) %>%
      summarise(min := min({{col_name}}),
                med := median({{col_name}}),
                max := max({{col_name}})) %>%
      ggplot(aes(sparseLevel, med, color = method, group = interaction(sparseLevel, method))) +
      # geom_boxplot(lwd = 0.5) +
      geom_point(size = 0.7, position = position_dodge(width = 1)) +
      geom_errorbar(aes(ymin = min, ymax = max), position = position_dodge(width = 1)) +
      facet_wrap(~ N, nrow = 1) +
      scale_color_manual(values = COLORVALUES) + 
      scale_y_continuous(breaks = seq(0,1,0.5), limits = c(0,1)) + 
      theme_classic() +
      guides(color = 'none') + 
      xlab("Relative sparsity") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    return(p)
  }
  
  
  mySummarise(F1_site) + ylab("F1 score")
  ggsave(here(plot_dir, paste0(experiment, "_chr_f1_siteLevel.png")), width = 4.5, height = 3)
  mySummarise(1-fdr_site) + ylab("Precision")
  ggsave(here(plot_dir, paste0(experiment, "_chr_precision_siteLevel.png")), width = 4.5, height = 3)
  mySummarise(power_site) + ylab("Recall")
  ggsave(here(plot_dir, paste0(experiment, "_chr_recall_siteLevel.png")), width = 4.5, height = 3)
  
  mySummarise(F1_region) + ylab("F1 score")
  ggsave(here(plot_dir, paste0(experiment, "_chr_f1_regionLevel.png")), width = 4.5, height = 3)
  mySummarise(1-fdr_region) + ylab("Precision")
  ggsave(here(plot_dir, paste0(experiment, "_chr_precision_regionLevel.png")), width = 4.5, height = 3)
  mySummarise(power_region) + ylab("Recall")
  ggsave(here(plot_dir, paste0(experiment, "_chr_recall_regionLevel.png")), width = 4.5, height = 3)
  
}

plotDefaultSetting(experiment = 'real')
plotDefaultSetting(experiment = 'sim')


# ---- RRA (AUPRC) ----

computeRRA <- function(power, fdr, xlim = c(0, 0.2)) { # The Ratio of Relevant Areas (RRA) indicator
  
  idx <- which(fdr >= xlim[1] & fdr <= xlim[2])
  pfdr <- c(xlim[1], fdr[idx], xlim[2])
  ppower <- c(power[idx], power[idx[length(idx)]])
  pauc <- sum(ppower * diff(pfdr))
  rra <- pauc / (xlim[2] - xlim[1])
  
  return(rra)
} 



plotRRA <- function(experiment) {
  
  getRRA <- function(N) {
    
    computeRRAbyExpr <- function(N, NP, sparseLevel) {
      
      smr0 <- fread(here(read_dir, 
                         paste0("benchmark_", experiment, "_chr"), 
                         "summary", 
                         paste0(ifelse(experiment == "real", "pseudo", "sim"),
                                "Chr_fdr&power_siteLevel_",
                                N , "cells_", NP, "subpops_",
                                NV, "VMRs_sparseLevel", sparseLevel,
                                "_seed", seed, ".csv.gz")))
      
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
      mutate(N = factor(paste0('N cells = ', N), levels = c('N cells = 200', 'N cells = 500', 'N cells = 1000', 'N cells = 2000')),
             sparseLevel = factor(sparseLevel, labels = c('high', 'medium', 'low')),
             method = factor(sapply(method, methodName), levels = LEVELS))
    return(smr_rra)
  }
  
  smr_rra0 <- rbind(
    getRRA(N = 200),
    getRRA(N = 500),
    getRRA(N = 1000),
    getRRA(N = 2000)
  ) 
  
  mySummarise <- function(col_name) {
    p <- smr_rra0 %>%
      group_by(N, sparseLevel, method) %>%
      summarise(min := min({{col_name}}),
                med := median({{col_name}}),
                max := max({{col_name}})) %>%
      ggplot(aes(sparseLevel, med, color = method, group = interaction(sparseLevel, method))) +
      # geom_boxplot(lwd = 0.5) +
      geom_point(size = 0.7, position = position_dodge(width = 1)) +
      geom_errorbar(aes(ymin = min, ymax = max), position = position_dodge(width = 1)) +
      facet_wrap(~ N, nrow = 1) +
      scale_color_manual(values = COLORVALUES) + 
      scale_y_continuous(breaks = seq(0,1,0.5), limits = c(0,1)) + 
      theme_classic() + 
      guides(color = 'none') + 
      xlab("Relative sparsity") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    return(p)
  }
  
  
  mySummarise(rra_site) + ylab("Ratio of relative areas")
  ggsave(here(plot_dir, paste0(experiment, "_chr_rra_siteLevel.png")), width = 4.5, height = 3)
  mySummarise(rra_region) + ylab("Ratio of relative areas")
  ggsave(here(plot_dir, paste0(experiment, "_chr_rra_regionLevel.png")), width = 4.5, height = 3)
  
}

plotRRA(experiment = 'real')
plotRRA(experiment = 'sim')


# ---- Distribution of spiked-in VMRs in simulation studies ----
summarySpikedInVMRs <- function(N, NP, sparseLevel) {
  dir <- paste0("data/interim/sim_studies/benchmark_real_chr/modified_real/pseudoChr_",
                subtype, "_", chromosome, "_",
                N, "cells_", NP, "subpops_",
                NV, "VMRs_sparseLevel", sparseLevel,
                "_seed", seed)
  gr <- granges(loadHDF5SummarizedExperiment(dir))
  h5closeAll()
  idx_vmr_start <- which(gr$is_vml==TRUE & lag(gr$is_vml)==FALSE)
  idx_vmr_end <- which(gr$is_vml==FALSE & lag(gr$is_vml)==TRUE) - 1
  vmr_n_cpg <- idx_vmr_end - idx_vmr_start
  vmr_width <- start(gr)[idx_vmr_end] - start(gr)[idx_vmr_start]
  return(data.frame(N, NP, sparseLevel, vmr_n_cpg, vmr_width))
}

Ns <- c(200, 500, 1000, 2000)
sparseLevels <- 1:3
params.df <- expand_grid(Ns, NP, sparseLevels)
smr.df <- with(params.df,  
               map_dfr(
                 1:nrow(params.df), 
                 ~ summarySpikedInVMRs(Ns[.x], 2, sparseLevels[.x])
               ))

smr.df %>%  
  mutate(
    N = factor(paste0('N cells = ', N), levels = c('N cells = 200', 'N cells = 500', 'N cells = 1000', 'N cells = 2000')),
    sparseLevel = factor(sparseLevel, labels = c('high', 'medium', 'low'))
  ) %>% 
  ggplot(aes(sparseLevel, vmr_width)) +
  geom_boxplot() +
  facet_wrap(~ N, nrow = 1) +
  scale_y_continuous(labels = scales::comma, 
                     limits = c(0, 25000), 
                     name = "Width of spiked-in VMRs in synthetic dataset (bp)") +
  xlab("Relative sparsity") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
ggsave(here(plot_dir, paste0("boxplot_vmr_width_in_simulations.png")), width = 6, height = 5)



