source("code/SETPATHS.R")

plot_dir <- "plots/sim_studies/benchmark_real_chr/"
NV <- 2000

# ---- gather results ----
time.df <- expand.grid(Method = c('vseq', 'vseq_cr', 'scbs', 'scmet'),
                       N = c(200, 500, 1000, 2000),
                       NP = c(2,3,4,5,8,12,20),
                       sparseLevel = 1:3,
                       time = 0)
for (i in 1:nrow(time.df)) {
  method <- time.df$Method[i]
  if (method == 'vseq_cr') {
    time.df$time[i] <- fread(paste0("data/interim/sim_studies/benchmark_real_chr/vmrseq/output/modelTime_crs_",
                                    time.df$N[i], "cells_", time.df$NP[i], "subpops_",
                                    NV, "VMRs_sparseLevel", time.df$sparseLevel[i],
                                    "_alpha", alpha, ".txt"))$time
    
  } else if (method == 'vseq') {
    alpha <- 0.05
    time.df$time[i] <- fread(paste0("data/interim/sim_studies/benchmark_real_chr/vmrseq/output/modelTime_",
                                    time.df$N[i], "cells_", time.df$NP[i], "subpops_",
                                    NV, "VMRs_sparseLevel", time.df$sparseLevel[i],
                                    "_alpha", alpha, ".txt"))$time
  } else if (method == 'scbs') {
    seed <- 2022
    time.df$time[i] <- fread(paste0("data/interim/sim_studies/benchmark_real_chr/scbs/output/modelTime_",
                                    time.df$N[i], "cells_", time.df$NP[i], "subpops_",
                                    NV, "VMRs_sparseLevel", time.df$sparseLevel[i],
                                    "_seed", seed, ".txt"))$time
  } else if (method == 'scmet') {
    bp_size <- 20000; seed <- 2022
    time.df$time[i] <- fread(paste0("data/interim/sim_studies/benchmark_real_chr/scmet/output/modelTime_", 
                                    bp_size/1000, "kbWindow_", 
                                    time.df$N[i], "cells_", time.df$NP[i], "subpops_", 
                                    NV, "VMRs_sparseLevel", time.df$sparseLevel[i],
                                    "_seed", seed, ".txt"))$time
  }
}

# ---- plotting ----

methodName <- function(method) switch (as.character(method),
                                       'vseq' = 'vmrseq',
                                       'vseq_cr' = 'vmrseq CRs',
                                       'scbs' = 'scbs',
                                       # 'smwd' = 'Smallwood',
                                       'scmet' = 'scMET')
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = COLORS[1], "vmrseq CRs" = COLORS[2],
                 "scbs" = COLORS[3], 
                 # "Smallwood" = COLORS[4], 
                 "scMET" = COLORS[5])

time.df %>%
  mutate(sparseLevel = as.factor(sparseLevel),
         Method = map_chr(Method, methodName),
         sparseLevel = factor(sparseLevel, labels = c('Relative sparsity = high', 
                                                      'Relative sparsity = medium', 
                                                      'Relative sparsity = low')),
  ) %>%
  ggplot(aes(N, time, 
             color = Method, 
             group = Method)) +
  # geom_jitter(width = 50) +
  geom_point() +
  geom_line(aes(group = interaction(Method, NP)), linetype = 'dashed', alpha = 0.5) +
  facet_wrap(~ sparseLevel) + 
  scale_color_manual(values = COLORVALUES) +
  scale_x_continuous(breaks = c(200, 500, 1000, 2000), labels = scales::comma) +
  xlab('# cells') + ylab('Time (min)') + 
  theme_classic()
ggsave(here(plot_dir, "point_scalability.png"), height = 4, width = 9)
