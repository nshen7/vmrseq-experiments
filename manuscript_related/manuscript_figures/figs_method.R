source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)

plot_dir <- "manuscript_related/manuscript_figures/method"
if (!file.exists(plot_dir)) dir.create(plot_dir)

COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = COLORS[1], "vmrseq CRs" = COLORS[2],
                 "scbs" = COLORS[3], "Smallwood" = COLORS[4], "scMET" = COLORS[5])

# ---- Illustration of smoothing ----
set.seed(2021)
plot.df <- data.frame(
  pos = cumsum(sample(50:100, 80, replace = T)),
  var = c(rnorm(20, 0.02, 0.02), 
          rnorm(3, 0.05, 0.02), rnorm(9, 0.1, 0.02), rnorm(3, 0.05, 0.02), 
          rnorm(25, 0.02, 0.02),
          rnorm(7, 0.04, 0.01),
          rnorm(13, 0.02, 0.02))
) %>%
  mutate(smoothed_var = lowess(pos, var, f = 1/15)$y)
cut <- 0.031
# (ind <- which(plot.df$smoothed_var >= cut))
# #  22 23 24 25 26 27 28 29 30 31 32 33 34 60 61 62 64 65 66 67 76 77 78
# start1 <- plot.df$pos[ind[1]]
# end1   <- plot.df$pos[ind[which(ind == 34)]]
# start2 <- plot.df$pos[ind[which(ind == 60)]]
# end2   <- plot.df$pos[ind[which(ind == 67)]]
plot.df %>%
  ggplot(aes(pos, smoothed_var)) + 
  # geom_rect(xmin = start1-50, xmax = end1+50, ymin = -Inf, ymax = Inf, alpha = 0.4, fill = COLORVALUES[3]) + 
  # geom_rect(xmin = start2-50, xmax = end2+50, ymin = -Inf, ymax = Inf, alpha = 0.4, fill = COLORVALUES[3]) + 
  geom_hline(aes(yintercept = cut), color = COLORVALUES[1], linetype = 'dashed', size = 2) +
  geom_point(size = 2) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank())
ggsave(here(plot_dir, "smooth_candidate_region.png"), width = 11, height = 3)


# ---- Illustration of underlying states for two-grouping case ----

x <- c(rep(1,17))
y <- c(rep(0,4),rep(1,9),rep(0,4))
pos <- c(1:4,4:12,12:15)

ggplot() +
  geom_path(aes(pos, x), color = "darkgoldenrod3", size = 2, alpha = 0.7) +
  geom_point(aes(pos[-c(4,14)], x[-c(4,14)]), shape = 21, size = 5) +
  ylim(-0.2, 1.2) +
  theme_void() 
ggsave(here(plot_dir, "states_grouping1.png"), width = 5, height = 1)
ggplot() +
  geom_path(aes(pos, y), color = "darkolivegreen3", size = 2, alpha = 1) +
  geom_point(aes(pos[-c(1:4,14:17)], y[-c(1:4,14:17)]), shape = 21, size = 5) +
  ylim(-0.2, 1.2) +
  # scale_fill_manual(values = c('white', '#D37538')) + 
  guides(fill = 'none') + 
  theme_void() 
ggsave(here(plot_dir, "states_grouping2.png"), width = 5, height = 1.1)

# ---- Illustration of underlying states for one-grouping case ----

x <- c(rep(1,9))
y <- c(rep(1,5),rep(0,4))
pos <- c(1:5,5:8)

ggplot() +
  geom_path(aes(pos, y), color = "grey", size = 2) +
  geom_point(aes(pos[-5], y[-5]), shape = 21, size = 5.5) +
  ylim(-0.2, 1.2) +
  guides(fill = 'none') + 
  theme_void() 
ggsave(here(plot_dir, "states_grouping_null.png"), width = 2.5, height = 1.1)


# ---- Illustration of groupings ----

x <- c(rep(0,7),rep(0,8),rep(0,7))
z <- y <- c(rep(0,7),rep(1,8),rep(0,7))
pos <- c(1:7,7:14,14:20)

png(here(plot_dir, "subpop1.png"), width = 300, height = 300, res = 70)
plot(pos, x, type = "l", ylim = c(0,1.4), ylab = "", xlab = "",
     col = "grey14", lwd = 7, yaxt = "n", xlim = c(0,20), bty="n", xaxt = "n")
# points(pos, x, col = "red", pch = 16)
axis(2, at = c(0,1))
legend(x = -1, y = 1.4, legend = "(Cell proportion = 0.3)", bty="n", cex = 1.5)
dev.off()
png(here(plot_dir, "subpop2.png"), width = 300, height = 300)
plot(pos, y, type = "l", ylim = c(0,1.4), ylab = "", xlab = "",
     col = "grey14", lwd = 7, yaxt = "n", xlim = c(0,20), bty="n", xaxt = "n")
# points(pos, y, col = "green", pch = 16)
axis(2, at = c(0,1))
legend(x = -1, y = 1.4, legend = "(Cell proportion = 0.2)", bty="n", cex = 1.5)
dev.off()
png(here(plot_dir, "subpop3.png"), width = 300, height = 300)
plot(pos, z, type = "l", ylim = c(0,1.4), ylab = "", xlab = "",
     col = "grey14", lwd = 7, yaxt = "n", xlim = c(0,20), bty="n", xaxt = "n")
# points(pos, z, col = "blue", pch = 16)
axis(2, at = c(0,1))
legend(x = -1, y = 1.4, legend = "(Cell proportion = 0.5)", bty="n", cex = 1.5)
dev.off()

png(here(plot_dir, "grouping1.png"), width = 300, height = 300)
plot(pos, x, type = "l", ylim = c(0,1.4), ylab = "", xlab = "",
     col = "darkgoldenrod3", lwd = 7, yaxt = "n", xlim = c(0,20), bty="n", xaxt = "n")
# points(pos, 0.3*x, col = "darkred", pch = 16)
axis(2, at = c(0,1))
legend(x = -1, y = 1.4, legend = "(Cell proportion = 0.3)", bty="n", cex = 1.5)
dev.off()

png(here(plot_dir, "grouping2.png"), width = 300, height = 300)
plot(pos, y, type = "l", ylim = c(0,1.4), ylab = "", xlab = "",
     col = "darkolivegreen3", lwd = 7, yaxt = "n", xlim = c(0,20), bty="n", xaxt = "n")
# points(pos, y, col = "darkgreen", pch = 16)
axis(2, at = c(0,1))
legend(x = -1, y = 1.4, legend = "(Cell proportion = 0.7)", bty="n", cex = 1.5)
dev.off()

png(here(plot_dir, "all_cells.png"), width = 300, height = 300)
plot(pos, 0.3*x + 0.2*y + 0.5*z, type = "l", ylim = c(0,1.2), ylab = "", xlab = "",
     col = "grey14", lwd = 7, yaxt = "n", xlim = c(0,20), bty="n", xaxt = "n")
# points(pos, 0.3*x + 0.2*y + 0.5*z, col = "black", pch = 16)
segments(0, 0.7, 7, 0.7, lty = 2, col = "grey")
axis(2, at = c(0,0.7,1))
dev.off()

png(here(plot_dir, "all_cells_null.png"), width = 300, height = 300)
plot(pos, 0*x, type = "l", ylim = c(0,1.2), ylab = "", xlab = "",
     col = "grey14", lwd = 7, yaxt = "n", xlim = c(0,20), bty="n", xaxt = "n")
axis(2, at = c(0,1))
dev.off()


# ---- Number of VMLs in hom study ----
read_dir <- 'data/interim/case_studies/luo2017mice_subset_hom/result_summary/'
sites.gr <- readRDS(paste0(read_dir, 'cpg_sites.rds'))
res_region <- list(
  vseq = loadHDF5SummarizedExperiment(paste0(read_dir, 'vmrseq_regionSummary_vmrs')),
  vseq_cr = loadHDF5SummarizedExperiment(paste0(read_dir, 'vmrseq_regionSummary_crs')),
  scbs = loadHDF5SummarizedExperiment(paste0(read_dir, 'scbs_regionSummary_vmrs')),
  smwd = loadHDF5SummarizedExperiment(paste0(read_dir, 'smallwood_regionSummary_vmrs'))
)
res.gr <- map(res_region, granges)

methods <- c('vmrseq', 'vmrseq CRs', 'scbs', 'Smallwood')
methods <- factor(methods, levels = methods)
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = COLORS[1], "vmrseq CRs" = COLORS[2],
                 "scbs" = COLORS[3], "Smallwood" = COLORS[4], "scMET" = COLORS[5])

# Percentage of sites in detected regions
n_sites <- sapply(res_region, function(se) findOverlaps(GenomicRanges::reduce(granges(se)), sites.gr) %>% length())
pct <- n_sites / length(sites.gr)
ggplot() + 
  geom_bar(aes(x = methods, y = pct, color = methods, fill = methods), stat = 'identity') + 
  scale_y_continuous(labels = percent, limits = c(0, 0.118), name = 'Percentage of CpGs in detected VMRs') +
  scale_color_manual(values = COLORVALUES) + 
  scale_fill_manual(values = COLORVALUES) + 
  xlab('Methods') + 
  theme_classic() +
  theme(legend.position = 'none')
ggsave(here(plot_dir, 'barplot_pctSites_vs_methods.png'), height = 3, width = 3)  

