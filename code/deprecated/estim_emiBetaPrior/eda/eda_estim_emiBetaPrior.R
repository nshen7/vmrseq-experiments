library(tidyverse)
setwd(here::here())
source("code/estim_emiBetaPrior/helper_functions.R")

dens_liu <- fread("data/interim/estim_emiBetaPrior/estimated_density_Liu2021_8subtypes.txt") %>%
  select(-c(Sample, GEO_accession, MedFinalReads)) %>%
  mutate(DataSource = "Liu2021")
dens_luo <- fread("data/interim/estim_emiBetaPrior/estimated_density_Luo2017_mice_6subtypes.txt") %>%
  select(-MedCovPct) %>%
  mutate(CellClass = replace(CellClass, CellClass=="Excitatory", "Exc")) %>%
  mutate(CellClass = replace(CellClass, CellClass=="Inhibitory", "Inh")) %>%
  mutate(DataSource = "Luo2017")

all(colnames(dens_liu)==colnames(dens_luo)) # check if the colnames are matched
dt_dens <- rbind(dens_liu, dens_luo)


# ==== EDA on density distribution ====
dt_dens %>% 
  ggplot(aes(x = Density.x, y = CumDens, color = CellClass, group = SubType)) + 
  geom_path(aes(linetype = DataSource)) +
  scale_color_brewer(palette = "Set2") +
  xlab("Value") +
  ylab("CDF of fractional methylation from subtypes") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_Liu2021&Luo2017.png", width = 7, height = 5)

# cdf of MF in lowly methylated sites
dt_dens %>% 
  filter(Density.x < 0.5) %>%
  group_by(CellClass, SubType, N_cell) %>%
  mutate(Density.y = Density.y / sum(Density.y)) %>%
  summarise(meanMF = sum(Density.y*Density.x)) %>%
  ggplot(aes(N_cell, meanMF)) +
  geom_point()
dt_dens %>% 
  filter(Density.x < 0.5) %>%
  group_by(SubType, N_cell) %>%
  mutate(Density.y = Density.y / sum(Density.y)) %>%
  mutate(CumDens = cumsum(Density.y)) %>%
  ggplot(aes(x = Density.x, y = CumDens, color = CellClass, group = SubType)) + 
  geom_path(aes(linetype = DataSource)) +
  scale_color_brewer(palette = "Set2") +
  xlab("Value") +
  ylab("CDF of fractional methylation from subtypes") +
  xlim(0,0.5) + ylim(0,1) +  
  ggtitle("Low methylated sites (MF < 0.5)") + 
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_coloredByCellType_unmethPop_Liu2021&Luo2017.png", width = 7, height = 5)
dt_dens %>% 
  filter(Density.x < 0.5) %>%
  group_by(SubType, N_cell) %>%
  mutate(Density.y = Density.y / sum(Density.y)) %>%
  mutate(CumDens = cumsum(Density.y)) %>%
  ggplot(aes(x = Density.x, y = CumDens, color = N_cell, group = SubType)) + 
  geom_path(aes(linetype = DataSource)) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  xlab("Value") +
  ylab("CDF of fractional methylation from subtypes") +
  xlim(0,0.5) + ylim(0,1) +  
  ggtitle("Low methylated sites (MF < 0.5)") + 
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_coloredByNcell_unmethPop_Liu2021&Luo2017.png", width = 7, height = 5)

# cdf of MF in highly methylated sites
dt_dens %>% 
  filter(Density.x > 0.5) %>%
  group_by(CellClass, SubType, N_cell) %>%
  mutate(Density.y = Density.y / sum(Density.y)) %>%
  summarise(meanMF = sum(Density.y*Density.x)) %>%
  ggplot(aes(N_cell, meanMF)) +
  geom_point()
dt_dens %>% 
  filter(Density.x > 0.5) %>%
  group_by(SubType, N_cell) %>%
  mutate(Density.y = Density.y / sum(Density.y)) %>%
  mutate(CumDens = cumsum(Density.y)) %>%
  ggplot(aes(x = Density.x, y = CumDens, color = CellClass, group = SubType)) + 
  geom_path(aes(linetype = DataSource)) +
  scale_color_brewer(palette = "Set2") +
  xlab("Value") +
  ylab("CDF of fractional methylation from subtypes") +
  xlim(0.5, 1) + ylim(0,1) +
  ggtitle("High methylated sites (MF > 0.5)") + 
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_coloredByCellType_methPop_Liu2021&Luo2017.png", width = 7, height = 5)
dt_dens %>% 
  filter(Density.x > 0.5) %>%
  group_by(SubType, N_cell) %>%
  mutate(Density.y = Density.y / sum(Density.y)) %>%
  mutate(CumDens = cumsum(Density.y)) %>%
  ggplot(aes(x = Density.x, y = CumDens, color = N_cell, group = SubType)) + 
  geom_path(aes(linetype = DataSource)) +
  scale_color_gradientn(colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous")) + 
  xlab("Value") +
  ylab("CDF of fractional methylation from subtypes") +
  xlim(0.5, 1) + ylim(0,1) +
  ggtitle("High methylated sites (MF > 0.5)") + 
  theme_bw() 
ggsave("plots/estim_emiBetaPrior/eda/cdf_MFdistribution_coloredByNcell_methPop_Liu2021&Luo2017.png", width = 7, height = 5)

## mean MF vs. N cells
dt_dens %>%
  filter(Density.x < 0.5) %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  mutate(Density.y = Density.y/sum(Density.y)) %>%
  summarise(MeanMF = sum(Density.x * Density.y)) %>%
  ggplot(aes(N_cell, MeanMF)) +
  geom_point(aes(color = DataSource, shape = CellClass)) +
  geom_smooth(method = "lm")
dt_dens %>%
  filter(Density.x > 0.5) %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  mutate(Density.y = Density.y/sum(Density.y)) %>%
  summarise(MeanMF = sum(Density.x * Density.y)) %>%
  ggplot(aes(N_cell, MeanMF)) +
  geom_point(aes(color = DataSource, shape = CellClass)) +
  geom_smooth(method = "lm")


## MF density on 0 / 1 vs. N cells
dt_dens %>%
  filter(Density.x < 0.5) %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  mutate(Density.y = Density.y/sum(Density.y)) %>%
  filter(Density.x == 0) %>%
  ggplot(aes(N_cell, Density.y)) +
  geom_point(aes(color = DataSource, shape = CellClass)) +
  geom_smooth(method = "lm")
dt_dens %>%
  filter(Density.x > 0.5) %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  mutate(Density.y = Density.y/sum(Density.y)) %>%
  filter(Density.x == 1) %>%
  ggplot(aes(N_cell, Density.y)) +
  geom_point(aes(color = DataSource, shape = CellClass)) +
  geom_smooth(method = "lm")


## fitted parameters for unmeth subpop vs. N cells
dt_dens %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  summarise(param1_u = mean(param1_u)) %>%
  ggplot(aes(N_cell, param1_u)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  ylab("beta_1") + xlab("N cells") + ylim(0,20)
ggsave("plots/estim_emiBetaPrior/eda/point_fittedPrior_unmethPop_param1_vs_Ncell_Liu2021&Luo2017.png", width = 5, height = 4)
dt_dens %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  summarise(param2_u = mean(param2_u)) %>%
  ggplot(aes(N_cell, param2_u)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) + 
  ylab("beta_2") + xlab("N cells") + ylim(1000,2250)
ggsave("plots/estim_emiBetaPrior/eda/point_fittedPrior_unmethPop_param2_vs_Ncell_Liu2021&Luo2017.png", width = 5, height = 4)
dt_dens %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  summarise(w1_u = mean(w1_u)) %>%
  ggplot(aes(N_cell, w1_u)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  ylab("weight") + xlab("N cells") + ylim(0.85, 1)
ggsave("plots/estim_emiBetaPrior/eda/point_fittedPrior_unmethPop_w1_vs_Ncell_Liu2021&Luo2017.png", width = 5, height = 4)
dt_dens %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  summarise(l1_diff_u = mean(l1_diff_u)) %>%
  ggplot(aes(N_cell, l1_diff_u)) +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) +
  ylab("L1 difference") + xlab("N cells") + ylim(0.25,0.5)
ggsave("plots/estim_emiBetaPrior/eda/point_fittedPrior_unmethPop_L1Diff_vs_Ncell_Liu2021&Luo2017.png", width = 5, height = 4)

## fitted parameters for meth subpop vs. N cells
dt_dens %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  summarise(param1_m = mean(param1_m)) %>%
  ggplot(aes(N_cell, param1_m)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2)  + 
  ylab("alpha_1") + xlab("N cells") + ylim(0,20)
ggsave("plots/estim_emiBetaPrior/eda/point_fittedPrior_methPop_param1_vs_Ncell_Liu2021&Luo2017.png", width = 5, height = 4)
dt_dens %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  summarise(param2_m = mean(param2_m)) %>%
  ggplot(aes(N_cell, param2_m)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) + 
  ylab("alpha_2") + xlab("N cells") + ylim(1000,2250)
ggsave("plots/estim_emiBetaPrior/eda/point_fittedPrior_methPop_param2_vs_Ncell_Liu2021&Luo2017.png", width = 5, height = 4)
dt_dens %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  summarise(w1_m = mean(w1_m)) %>%
  ggplot(aes(N_cell, w1_m)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) + 
  ylab("weight") + xlab("N cells") + ylim(0.85, 1)
ggsave("plots/estim_emiBetaPrior/eda/point_fittedPrior_methPop_w1_vs_Ncell_Liu2021&Luo2017.png", width = 5, height = 4)
dt_dens %>%
  group_by(SubType, N_cell, CellClass, DataSource) %>%
  summarise(l1_diff_m = mean(l1_diff_m)) %>%
  ggplot(aes(N_cell, l1_diff_m)) +
  geom_point(aes(color = DataSource, shape = CellClass), size = 2) + 
  ylab("L1 difference") + xlab("N cells") + ylim(0.25,0.5)
ggsave("plots/estim_emiBetaPrior/eda/point_fittedPrior_methPop_L1Diff_vs_Ncell_Liu2021&Luo2017.png", width = 5, height = 4)

# ==== Effect of cell types on prior estimation (Investigate two subtypes that seems far apart: DG dg-all and PAL-Inh Meis2) ====
x <- seq(from = 0, to = 1, by = 1/n_bins)

best_params_u_1 <- dt_dens %>%
  filter(SubType == "DG dg-all") %>%
  select(1,2) %>%
  calBestParamFromDensity(population = "u", mc.cores = 16)
# Fitted mixture distribution density for unmethylated population: 0.61 * dbeta(x, 1, 7) + 0.39 * dbeta(x, 1, 210) 
u1 <- 0.61 * pbeta(x, 1, 7) + 0.39 * pbeta(x, 1, 210)

best_params_m_1 <- dt_dens %>%
  filter(SubType == "DG dg-all") %>%
  select(1,2) %>%
  calBestParamFromDensity(population = "m", mc.cores = 16)
# Fitted mixture distribution density for methylated population: 0.99 * dbeta(x, 13, 1) + 0.01 * dbeta(x, 430, 1) 
m1 <- 0.99 * pbeta(x, 13, 1) + 0.01 * pbeta(x, 430, 1)

best_params_u_2 <- dt_dens %>%
  filter(SubType == "PAL-Inh Meis2") %>%
  select(1,2) %>%
  calBestParamFromDensity(population = "u", mc.cores = 16)
# Fitted mixture distribution density for unmethylated population: 0.88 * dbeta(x, 1, 5) + 0.12 * dbeta(x, 1, 1760)
u2 <- 0.88 * pbeta(x, 1, 5) + 0.12 * pbeta(x, 1, 1760)

best_params_m_2 <- dt_dens %>%
  filter(SubType == "PAL-Inh Meis2") %>%
  select(1,2) %>%
  calBestParamFromDensity(population = "m", mc.cores = 16)
# Fitted mixture distribution density for methylated population: 0.96 * dbeta(x, 12, 1) + 0.04 * dbeta(x, 1600, 1) 
m2 <- 0.96 * pbeta(x, 12, 1) + 0.04 * pbeta(x, 1600, 1) 


best_params_u_3 <- dt_dens %>%
  filter(SubType == "CA1 Chrm3") %>%
  select(1,2) %>%
  calBestParamFromDensity(population = "u", mc.cores = 16)
# Fitted mixture distribution density for unmethylated population: 0.93 * dbeta(x, 1, 6) + 0.07 * dbeta(x, 1, 1860) 
u3 <- 0.93 * pbeta(x, 1, 6) + 0.07 * pbeta(x, 1, 1860) 

best_params_m_3 <- dt_dens %>%
  filter(SubType == "CA1 Chrm3") %>%
  select(1,2) %>%
  calBestParamFromDensity(population = "m", mc.cores = 16)
# Fitted mixture distribution density for methylated population: 0.98 * dbeta(x, 12, 1) + 0.02 * dbeta(x, 1540, 1) 
m3 <- 0.98 * pbeta(x, 12, 1) + 0.02 * pbeta(x, 1540, 1) 

1-m1[n_bins]; 1-m2[n_bins]; 1-m3[n_bins]

### Compare DG dg-all &  PAL-Inh Meis2
data.frame(u1, u2) %>%
  ggplot(aes(u1, u2)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() +
  ggtitle("QQplot of fitted prior for unmethylated population from 2 subtypes") +
  xlab("Subtype DG dg-all") + ylab("Subtype PAL-Inh Meis2")
ggsave("plots/estim_emiBetaPrior/eda/qqplot_fittedUnmethPrior_subtypes_DG_dg-all_vs_PAL-Inh_Meis2.png", width = 6, height = 5)

data.frame(m1, m2) %>%
  ggplot(aes(m1, m2)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() +
  ggtitle("QQplot of fitted prior for methylated population from 2 subtypes") +
  xlab("Subtype DG dg-all") + ylab("Subtype PAL-Inh Meis2")
ggsave("plots/estim_emiBetaPrior/eda/qqplot_fittedMethPrior_subtypes_DG_dg-all_vs_PAL-Inh_Meis2.png", width = 6, height = 5)

dt_dens %>%
  filter(SubType %in% c("DG dg-all", "PAL-Inh Meis2")) %>%
  select(Density.x, Density.y, SubType) %>%
  pivot_wider(names_from = SubType, values_from = Density.y) %>%
  mutate(DG_cumDens = cumsum(`DG dg-all`), PAL_cumDens = cumsum(`PAL-Inh Meis2`)) %>%
  ggplot(aes(DG_cumDens, PAL_cumDens)) + 
  geom_point() +
  geom_point(x = 0, y = 0) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() +
  ggtitle("QQplot of across-cell MF from 2 subtypes") +
  xlab("Subtype DG dg-all") + ylab("Subtype PAL-Inh Meis2")
ggsave("plots/estim_emiBetaPrior/eda/qqplot_MF_subtypes_DG_dg-all_vs_PAL-Inh_Meis2.png", width = 6, height = 5)


### compare CA1 Chrm3 and PAL-Inh Meis2
data.frame(u3, u2) %>%
  ggplot(aes(u3, u2)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() +
  ggtitle("QQplot of fitted prior for unmethylated population from 2 subtypes") +
  xlab("Subtype CA1 Chrm3") + ylab("Subtype PAL-Inh Meis2")
ggsave("plots/estim_emiBetaPrior/eda/qqplot_fittedUnmethPrior_subtypes_CA1_Chrm3_vs_PAL-Inh_Meis2.png", width = 6, height = 5)

data.frame(m3, m2) %>%
  ggplot(aes(m3, m2)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() +
  ggtitle("QQplot of fitted prior for methylated population from 2 subtypes") +
  xlab("Subtype CA1 Chrm3") + ylab("Subtype PAL-Inh Meis2")
ggsave("plots/estim_emiBetaPrior/eda/qqplot_fittedMethPrior_subtypes_CA1_Chrm3_vs_PAL-Inh_Meis2.png", width = 6, height = 5)

dt_dens %>%
  filter(SubType %in% c("CA1 Chrm3", "PAL-Inh Meis2")) %>%
  select(Density.x, Density.y, SubType) %>%
  pivot_wider(names_from = SubType, values_from = Density.y) %>%
  mutate(CA1_cumDens = cumsum(`CA1 Chrm3`), PAL_cumDens = cumsum(`PAL-Inh Meis2`)) %>%
  ggplot(aes(CA1_cumDens, PAL_cumDens)) + 
  geom_point() +
  geom_point(x = 0, y = 0) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0,1.01) + ylim(0,1.01) +
  theme_bw() +
  ggtitle("QQplot of across-cell MF from 2 subtypes") +
  xlab("Subtype CA1 Chrm3") + ylab("Subtype PAL-Inh Meis2")
ggsave("plots/estim_emiBetaPrior/eda/qqplot_MF_subtypes_CA1_Chrm3_vs_PAL-Inh_Meis2.png", width = 6, height = 5)


# ==== Compare priors estimated in DXM and estimated by ourselves before comprehensive exam ====
x <- seq(from = 0, to = 1, by = 0.005)
m0 <- 0.1 * pbeta(x, 200, 1) + 0.9 * pbeta(x, 6, 1) # from DXM paper
u0 <- 0.61 * pbeta(x, 1, 100) + 0.39 * pbeta(x, 1, 3) # from DXM paper
m00 <- 0.08 * pbeta(x, 950, 1) + 0.92 * pbeta(x, 9, 1) # fitted in pre-comp experiments
u00 <- 0.13 * pbeta(x, 1, 100) + 0.87 * pbeta(x, 1, 3) # fitted in pre-comp experiments

data.frame(u0, u00) %>%
  ggplot(aes(u0, u00)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() +
  ggtitle("QQplot of fitted prior for unmeth pop from DXM and pre-comp expr") +
  xlab("From DXM") + ylab("From pre-comp experiments")
ggsave("plots/estim_emiBetaPrior/eda/qqplot_fittedUnmethPrior_DXMPrior_vs_preCompPrior.png", width = 6, height = 5)

data.frame(m0, m00) %>%
  ggplot(aes(m0, m00)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlim(0,1) + ylim(0,1) +
  theme_bw() +
  ggtitle("QQplot of fitted prior for meth pop from DXM and pre-comp expr") +
  xlab("From DXM") + ylab("From pre-comp experiments")
ggsave("plots/estim_emiBetaPrior/eda/qqplot_fittedMethPrior_DXMPrior_vs_preCompPrior.png", width = 6, height = 5)


  
  

  
  
  
  
  
  
  
  

