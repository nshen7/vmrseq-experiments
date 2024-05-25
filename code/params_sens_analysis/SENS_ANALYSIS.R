source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
plot_dir <- "code/params_sens_analysis"
  
# ---- beta_prior_m ----
fit_m_all <- fread("../vmrseq-package/vmrseq/data-raw/fitted_params_beta_prior_m.csv")
fit_m_liu <- fread("code/params_sens_analysis/fitted_params_beta_prior_m.csv")

ggplot() +
  geom_point(data = fit_m_liu, aes(log(med_cov), mu_fit), color = "orange") +
  geom_point(data = fit_m_all, aes(log(med_cov), mu_fit), color = "blue") +
  xlab("Log median across-cell coverage") + ylab("BB estimates of mu from all data and Liu only") +
  ylim(0.82, 0.96) +
  theme_classic()
ggsave(here(plot_dir, "point_methPop_paramMu_vs_medCov_allData&liuOnly.png"), width = 4, height = 4)
sub_meth_smr %>% ggplot() +
  geom_point(data = fit_m_liu, aes(log(med_cov), sigma_fit), color = "orange") +
  geom_point(data = fit_m_all, aes(log(med_cov), sigma_fit), color = "blue") +
  xlab("Log median across-cell coverage") + ylab("BB estimates of sigma from all data and Liu only") +
  ylim(0.065, 0.12) +
  theme_classic()
ggsave(here(plot_dir, "point_methPop_paramSigma_vs_medCov_allData&liuOnly.png"), width = 4, height = 4)

# ---- beta_prior_u ----
fit_u_all <- fread("../vmrseq-package/vmrseq/data-raw/fitted_params_beta_prior_u.csv")
fit_u_liu <- fread("code/params_sens_analysis/fitted_params_beta_prior_u.csv")

ggplot() +
  geom_point(data = fit_u_liu, aes(log(med_cov), nu_fit), color = "orange") +
  geom_point(data = fit_u_all, aes(log(med_cov), nu_fit), color = "blue") +
  xlab("Log median across-cell coverage") + ylab("BB estimates of nu from all data and Liu only") +
  ylim(-0.001, 0.6) +
  theme_classic()
ggsave(here(plot_dir, "point_unmethPop_paramNu_vs_medCov_allData&liuOnly.png"), width = 4, height = 4)
ggplot() +
  geom_point(data = fit_u_liu, aes(log(med_cov), mu_fit), color = "orange") +
  geom_point(data = fit_u_all, aes(log(med_cov), mu_fit), color = "blue") +
  xlab("Log median across-cell coverage") + ylab("BB estimates of mu from all data and Liu only") +
  ylim(0.05, 0.27) +
  theme_classic()
ggsave(here(plot_dir, "point_unmethPop_paramMu_vs_medCov_allData&liuOnly.png"), width = 4, height = 4)
sub_meth_smr %>% ggplot() +
  geom_point(data = fit_u_liu, aes(log(med_cov), sigma_fit), color = "orange") +
  geom_point(data = fit_u_all, aes(log(med_cov), sigma_fit), color = "blue") +
  xlab("Log median across-cell coverage") + ylab("BB estimates of sigma from all data and Liu only") +
  ylim(-0.001, 0.4) +
  theme_classic()
ggsave(here(plot_dir, "point_unmethPop_paramSigma_vs_medCov_allData&liuOnly.png"), width = 4, height = 4)


## get legend 
fit_u <- rbind(
  data.frame(fit_u_all, train_data = "All data"),
  data.frame(fit_u_liu, train_data = "Liu only")
) 
p <- fit_u %>% 
  ggplot(aes(log(med_cov), sigma_fit, color = train_data)) +
  geom_point() + 
  scale_color_manual(values = c("All data" = "blue", "Liu et al only" = "orange"), name = "Training data") +
  theme_classic()
legend <- cowplot::get_legend(p)
png(here(plot_dir, "legend_sens_point.png"), width = 400, height = 300, res = 350)
grid::grid.newpage()
grid::grid.draw(legend)
dev.off()


# ---- transition prob ----
tp_all <- vmrseq:::tp0
tp_liu <- readRDS('code/params_sens_analysis/tp0.rds')

plot_tp_all <- data.frame(
    dist_bp = 1:tp_all@max_dist_bp,
    tp_all@transit_probs,
    train_data = "All data"
) 
plot_tp_liu <- data.frame(
  dist_bp = 1:tp_liu@max_dist_bp,
  tp_liu@transit_probs,
  train_data = "Liu only"
) 
plot_tp <- rbind(plot_tp_all, plot_tp_liu) %>%
  tidyr::pivot_longer(
    cols = -c('dist_bp', 'train_data'),
    names_to = c(".value", "type"),
    names_pattern = "(.*)_(.*)"
  )

type_labs <- c("P(0|0)","P(0|1)","P(1|0)","P(1|1)")
names(type_labs) <- c("00","01","10","11")

plot_tp %>%
  ggplot() +
  geom_path(aes(dist_bp, phat, color = train_data)) +
  theme_bw() + ylim(0, 1) +
  scale_color_manual(values = c("All data" = "blue", "Liu only" = "orange"), name = "Training data") +
  xlab("CpG-CpG Distance") + ylab("Transition Probability") +
  facet_wrap(~ type, labeller = labeller(type = type_labs))
ggsave(here(plot_dir, "transit_probs_allData&liuOnly.png"), width = 6, height = 5)
