.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
devtools::load_all("../vmrseq-package/vmrseq/")
library(tidyverse)
library(data.table)
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(scales)
library(pheatmap)

read_dir <- "data/interim/case_studies/luo2017mice_subset_het/"
write_dir <- "data/interim/case_studies/luo2017mice_subset_het/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_subset_het/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

md <- fread("data/interim/case_studies/luo2017mice_subset_het/metadata_luo2017mice_subset_het.csv")
rownames(md) <- md$sample
table(md$Neuron_type3)

# ==== read in results from various methods ====
sites.gr <- readRDS(paste0(read_dir, "result_summary/cpg_sites.rds"))
res_region <- list(
  vseq = loadHDF5SummarizedExperiment(paste0(read_dir, "result_summary/vmrseq_regionSummary_vmrs")),
  vseq_cr = loadHDF5SummarizedExperiment(paste0(read_dir, "result_summary/vmrseq_regionSummary_crs")),
  scbs = loadHDF5SummarizedExperiment(paste0(read_dir, "result_summary/scbs_regionSummary_vmrs")),
  smwd = loadHDF5SummarizedExperiment(paste0(read_dir, "result_summary/smallwood_regionSummary_vmrs")),
  scmet = loadHDF5SummarizedExperiment(paste0(read_dir, "result_summary/scmet_regionSummary_vmrs"))
)
res_site <- list(
  vseq = loadHDF5SummarizedExperiment(paste0(read_dir, "result_summary/vmrseq_siteSummary_sparseRep_vmrs")),
  vseq_cr = loadHDF5SummarizedExperiment(paste0(read_dir, "result_summary/vmrseq_siteSummary_sparseRep_crs")),
  scbs = loadHDF5SummarizedExperiment(paste0(read_dir, "result_summary/scbs_siteSummary_sparseRep_vmrs")),
  smwd = loadHDF5SummarizedExperiment(paste0(read_dir, "result_summary/smallwood_siteSummary_sparseRep_vmrs")),
  scmet = loadHDF5SummarizedExperiment(paste0(read_dir, "result_summary/scmet_siteSummary_sparseRep_vmrs"))
)

methods <- c("vmrseq", "CR in vmrseq", "scbs", "smallwood", "scMET")
methods <- factor(methods, levels = methods)


# ==== Color settings ====
COLORS <- RColorBrewer::brewer.pal(n = 6, name = "PuOr")[-3]
COLORVALUES <- c("vmrseq" = COLORS[1], "CR in vmrseq" = COLORS[2], 
                 "scbs" = COLORS[3], "smallwood" = COLORS[4], "scMET" = COLORS[5])


############################
##### Plot example VMR #####
############################

plotTargetRegion <- function(regions, 
                             extend = NULL,
                             nn = 100) {
  
  stopifnot("Multiple chromosomes!" = length(unique(regions$seqnames)) <= 1)
  
  seqname0 <- as.character(regions$seqnames[1])
  start0 <- min(regions$start)
  end0 <- max(regions$end)
  
  if (is.null(extend)) extend <- (end0 - start0) / 2
  
  # Target region
  entire_region.gr <- GRanges(
    seqnames = seqname0, 
    ranges = IRanges(start = start0-extend, end = end0+extend)
  )
  gr <- granges(res_region$vseq); vseq_olap <- subset(gr, countOverlaps(gr, entire_region.gr)>0) %>% as.data.frame()
  gr <- granges(res_region$vseq_cr); vseq_cr_olap <- subset(gr, countOverlaps(gr, entire_region.gr)>0) %>% as.data.frame()
  gr <- granges(res_region$scbs); scbs_olap <- subset(gr, countOverlaps(gr, entire_region.gr)>0) %>% as.data.frame()
  gr <- granges(res_region$smwd); smwd_olap <- subset(gr, countOverlaps(gr, entire_region.gr)>0) %>% as.data.frame()
  gr <- granges(res_region$scmet); scmet_olap <- subset(gr, countOverlaps(gr, entire_region.gr)>0) %>% as.data.frame()

  # Load in site-level methylation
  all_sites.se <- loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq/input/", seqname0))
  
  # Extract site-level info of target region
  idx <- queryHits(findOverlaps(all_sites.se, entire_region.gr))
  target_sites.se <- all_sites.se[idx]
  
  # Summarize mean methyl level across cells in each subtype
  for (subtype in unique(md$Neuron_type3)) {
    cell_idx <- which(md$Neuron_type3==subtype)
    total_meth <- rowSums(assays(target_sites.se)$M_mat[, cell_idx]) %>% round 
    total_count <- rowSums(assays(target_sites.se)$M_mat[, cell_idx] > 0)
    df <- data.frame(total_meth / total_count, total_count)
    colnames(df) <- paste0(c("mf_", "cov_"), subtype)
    rowData(target_sites.se) <- cbind(rowData(target_sites.se), df)
  }
  
  smr <- granges(target_sites.se) %>%
    as.data.frame() %>%
    pivot_longer(cols = starts_with("mf_") | starts_with("cov_"),
                 names_to = c("value_type", "subtype"),
                 names_sep = "_") %>%
    pivot_wider(id_cols = !c(meth, total),
                names_from = value_type, 
                values_from = value) 
  p <- ggplot() +
    stat_smooth(data = smr, aes(start, mf, color = subtype, weight = cov), method = "loess", span = nn/nrow(target_sites.se)) +
    geom_jitter(data = smr, aes(start, mf, color = subtype, size = cov), alpha = 0.5, height = 0.02) + 
    scale_size(range = c(1, 4)) +
    scale_color_brewer(palette = "Set2") +
    scale_fill_manual(name = "Methods", values = COLORVALUES)  + 
    coord_cartesian(ylim = c(-0.1, 1), xlim = c(start0-extend, end0+extend)) +
    theme_classic() +
    ggtitle(paste0(seqname0, ": ", start0, " - ", end0))
  if (nrow(vseq_olap) > 0) p <- p + geom_rect(data = vseq_olap, aes(xmin = start, xmax = end, fill = "vmrseq"), ymin = -0.04, ymax = -0.02) 
  if (nrow(vseq_cr_olap) > 0) p <- p +  geom_rect(data = vseq_cr_olap, aes(xmin = start, xmax = end, fill = "CR in vmrseq"), ymin = -0.06, ymax = -0.04) 
  if (nrow(scbs_olap) > 0) p <- p +  geom_rect(data = scbs_olap, aes(xmin = start, xmax = end, fill = "scbs"), ymin = -0.08, ymax = -0.06) 
  if (nrow(smwd_olap) > 0) p <- p +  geom_rect(data = smwd_olap, aes(xmin = start, xmax = end, fill = "smallwood"), ymin = -0.1, ymax = -0.08) 
  if (nrow(scmet_olap) > 0) p <- p + geom_rect(data = scmet_olap, aes(xmin = max(start0, start), xmax = min(end0, end), fill = "scMET"), ymin = -0.12, ymax = -0.1) 
  p
  return(p)
}


top_vmrs <- granges(res_region$vseq) %>%
  as.data.frame() %>% 
  arrange(desc(loglik_diff)) 
View(top_vmrs)


# # ==== plot example regions ====
# (r1 <- top_vmrs[1:11,])
# #    seqnames    start      end width strand index num_cpg start_ind end_ind        pi loglik_diff
# # 1      chr1 56905863 56908566  2704      *   759      16    351601  351616 0.4966625    884.6262
# # 2      chr1 56909429 56911812  2384      *   760      18    351621  351638 0.4966625    884.6262
# # 3      chr1 56912613 56926246 13634      *   761      94    351640  351733 0.4966625    884.6262
# # 4      chr1 56926833 56927822   990      *   762      10    351742  351751 0.4966625    884.6262
# # 5      chr1 56928499 56929357   859      *   763      15    351754  351768 0.4966625    884.6262
# # 6      chr1 56930883 56932392  1510      *   764      14    351775  351788 0.4966625    884.6262
# # 7      chr1 56933370 56934708  1339      *   765      21    351795  351815 0.4966625    884.6262
# # 8      chr1 56935130 56940267  5138      *   766      38    351820  351857 0.4966625    884.6262
# # 9      chr1 56942184 56942805   622      *   767       6    351863  351868 0.4966625    884.6262
# # 10     chr1 56943321 56946752  3432      *   768      41    351871  351911 0.4966625    884.6262
# # 11     chr1 56947813 56968847 21035      *   769     197    351920  352116 0.4966625    884.6262
# r1$start - lag(r1$end)
# plotTargetRegion(regions = r1)
# ggsave(paste0(plot_dir, "example_region_1.png"), width = 12, height = 6)
# 
# (r2 <- top_vmrs[96:97,])
# #    seqnames    start      end width strand index num_cpg start_ind end_ind        pi loglik_diff
# # 96    chr11 12025104 12027096  1993      *   125     163     76132   76294 0.3028376    234.1303
# # 97    chr11 12027921 12028152   232      *   126       9     76296   76304 0.3028376    234.1303
# r2$start - lag(r2$end)
# plotTargetRegion(regions = r2)
# ggsave(paste0(plot_dir, "example_region_2.png"), width = 12, height = 6)
# 
# (r3 <- top_vmrs[290:293,])
# #     seqnames    start      end width strand index num_cpg start_ind end_ind        pi loglik_diff
# # 290    chr17 51737704 51738049   346      *   852       6    459260  459265 0.6089023    129.3271
# # 291    chr17 51739428 51740310   883      *   853      23    459276  459298 0.6089023    129.3271
# # 292    chr17 51740637 51741860  1224      *   854      17    459306  459322 0.6089023    129.3271
# # 293    chr17 51742255 51743494  1240      *   855      26    459329  459354 0.6089023    129.3271
# r3$start - lag(r3$end)
# plotTargetRegion(regions = r3)
# ggsave(paste0(plot_dir, "example_region_3.png"), width = 12, height = 6)
# 
# (r4 <- top_vmrs[2831:2832,])
# #      seqnames     start       end width strand index num_cpg start_ind end_ind        pi loglik_diff
# # 2831     chr3 149910168 149910713   546      *  1921      10   1048928 1048937 0.4171892    40.63232
# # 2832     chr3 149910990 149911733   744      *  1922      17   1048945 1048961 0.4171892    40.63232
# r4$start - lag(r4$end)
# plotTargetRegion(regions = r4)
# ggsave(paste0(plot_dir, "example_region_4.png"), width = 12, height = 6)
# 
# (r5 <- top_vmrs[340:343,])
# #     seqnames    start      end width strand index num_cpg start_ind end_ind        pi loglik_diff
# # 340    chr13 42681330 42682571  1242      *   455      21    280954  280974 0.4711736    122.2384
# # 341    chr13 42683056 42686278  3223      *   456      27    280978  281004 0.4711736    122.2384
# # 342    chr13 42689411 42691273  1863      *   457      23    281023  281045 0.4711736    122.2384
# # 343    chr13 42692382 42692763   382      *   458       7    281055  281061 0.4711736    122.2384
# r5$start - lag(r5$end)
# plotTargetRegion(regions = r5)
# ggsave(paste0(plot_dir, "example_region_5.png"), width = 12, height = 6)
# 
# (r6 <- top_vmrs[5202:5203,])
# #      seqnames     start       end width strand index num_cpg start_ind end_ind        pi loglik_diff
# # 5202    chr11 116822879 116823021   143      *  1753       5   1051973 1051977 0.4676285     30.2759
# # 5203    chr11 116823849 116824261   413      *  1754      11   1051980 1051990 0.4676285     30.2759
# r6$start - lag(r6$end)
# plotTargetRegion(regions = r6)
# ggsave(paste0(plot_dir, "example_region_6.png"), width = 12, height = 6)
# 


# ==== plot same set of example regions containg VMRs ====
vmrs.gr <- granges(res_region$vseq)

r1.gr <- GRanges(
  seqnames = "chr1", 
  ranges = IRanges(start = 56905863, end = 56968847)
) 
plotTargetRegion(regions = as.data.frame(vmrs.gr[countOverlaps(vmrs.gr, r1.gr)>0]))
ggsave(paste0(plot_dir, "example_region_1.png"), width = 12, height = 6)

r2.gr <- GRanges(
  seqnames = "chr11", 
  ranges = IRanges(start = 12025104, end = 12028152)
) 
plotTargetRegion(regions = as.data.frame(vmrs.gr[countOverlaps(vmrs.gr, r2.gr)>0]))
ggsave(paste0(plot_dir, "example_region_2.png"), width = 12, height = 6)

r3.gr <- GRanges(
  seqnames = "chr17", 
  ranges = IRanges(start = 51737704, end = 51743494)
) 
plotTargetRegion(regions = as.data.frame(vmrs.gr[countOverlaps(vmrs.gr, r3.gr)>0]))
ggsave(paste0(plot_dir, "example_region_3.png"), width = 12, height = 6)

r4.gr <- GRanges(
  seqnames = "chr3", 
  ranges = IRanges(start = 149910168, end = 149911733)
) 
plotTargetRegion(regions = as.data.frame(vmrs.gr[countOverlaps(vmrs.gr, r4.gr)>0]))
ggsave(paste0(plot_dir, "example_region_4.png"), width = 12, height = 6)

r5.gr <- GRanges(
  seqnames = "chr13",
  ranges = IRanges(start = 42681330, end = 42692763)
) 
plotTargetRegion(regions = as.data.frame(vmrs.gr[countOverlaps(vmrs.gr, r5.gr)>0]))
ggsave(paste0(plot_dir, "example_region_5.png"), width = 12, height = 6)

r6.gr <- GRanges(
  seqnames = "chr11",
  ranges = IRanges(start = 116822879, end = 116824261)
) 
plotTargetRegion(regions = as.data.frame(vmrs.gr[countOverlaps(vmrs.gr, r6.gr)>0]))
ggsave(paste0(plot_dir, "example_region_6.png"), width = 12, height = 6)


r7.gr <- GRanges(
  seqnames = "chr16",
  ranges = IRanges(start = 49804106, end = 49804762)
) 
plotTargetRegion(regions = as.data.frame(vmrs.gr[countOverlaps(vmrs.gr, r7.gr)>0]), extend = 10000)
ggsave(paste0(plot_dir, "example_region_7.png"), width = 12, height = 6)


# ==== plot some CRs that does not contain VMRs ====
vmrs.gr <- granges(res_region$vseq)
crs.gr <- granges(res_region$vseq_cr)
crs_novmr.gr <- crs.gr[countOverlaps(crs.gr, vmrs.gr) == 0]

r.gr <- GRanges(
  seqnames = "chr4",
  ranges = IRanges(start = 113739235, end = 113747793)
) 
plotTargetRegion(regions = as.data.frame(r.gr), extend = 20000)

r.gr <- GRanges(
  seqnames = "chr7",
  ranges = IRanges(start = 42425702, end = 42436954)
) 
plotTargetRegion(regions = as.data.frame(r.gr), extend = 20000)

# ==== plot subtype marker genes ====

genes_md <- read.csv("data/metadata/metadata_luo2017/NIHMS893063-supplement-Table_S3_csv.csv", skip = 1)[1:22]
names(genes_md)[7:22] <- genes_md[1,7:22]
genes_md <- genes_md[-1, ] %>% select(c(1:6, which(names(genes_md) %in% c("mDL-2", "mPv", "mSst-1", "mSst-2"))))

idx <- 4
r.gr <- GRanges(
  seqnames = paste0("chr", genes_md$chr[idx]),
  ranges = IRanges(start = genes_md$start[idx], end = genes_md$end[idx])
) 
plotTargetRegion(regions = as.data.frame(r.gr), extend = 20000)

idx <- 7
r.gr <- GRanges(
  seqnames = paste0("chr", genes_md$chr[idx]),
  ranges = IRanges(start = genes_md$start[idx], end = genes_md$end[idx])
) 
plotTargetRegion(regions = as.data.frame(r.gr), extend = 20000)

idx <- 23
r.gr <- GRanges(
  seqnames = paste0("chr", genes_md$chr[idx]),
  ranges = IRanges(start = genes_md$start[idx], end = genes_md$end[idx])
) 
plotTargetRegion(regions = as.data.frame(r.gr), extend = 20000)

idx <- 26
r.gr <- GRanges(
  seqnames = paste0("chr", genes_md$chr[idx]),
  ranges = IRanges(start = genes_md$start[idx], end = genes_md$end[idx])
) 
plotTargetRegion(regions = as.data.frame(r.gr), extend = 20000)

                  