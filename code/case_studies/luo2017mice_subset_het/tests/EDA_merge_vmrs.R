read_dir <- "data/interim/case_studies/luo2017mice_subset_het/"
dirs_vseq <- paste0(read_dir, "vmrseq/output/chr", 1:19, ".rds")
dirs_vseq <- paste0(read_dir, "vmrseq/output/chr", 1:19, "_1kbWindow.rds")
gr <- do.call(c, lapply(dirs_vseq, function(dir) readRDS(dir)$gr))
vmrs.gr <- do.call(c, lapply(dirs_vseq, function(dir) readRDS(dir)$vmr.ranges))
crs.gr <- do.call(c, lapply(dirs_vseq, function(dir) readRDS(dir)$cr.ranges))

gr_cr <- gr %>% subset(!is.na(cr_index))
gr_cr$cr_index <- paste0(seqnames(gr_cr), "_", gr_cr$cr_index)
gr_cr$vmr_index <- paste0(seqnames(gr_cr), "_", gr_cr$vmr_index)
View(data.frame(gr_cr))

smr_vmr <- gr_cr %>% as_tibble() %>%
  filter(!is.na(cr_index)) %>%
  group_by(vmr_index) %>%
  summarise(n = n(), cr_index = max(cr_index))
n_vmrs_in_cr <- table(smr_vmr$cr_index)
table(n_vmrs_in_cr)
sum(n_vmrs_in_cr > 1) / length(n_vmrs_in_cr)


gr_cr_multi_vmr <- gr_cr %>% subset(cr_index %in% names(n_vmrs_in_cr)[n_vmrs_in_cr>1])
View(data.frame(gr_cr_multi_vmr))

vmrs.gr %>% as.data.frame() %>% arrange(desc(loglik_diff)) %>% View
quantile(vmrs.gr$num_cpg)
entire_region.gr <- GRanges(
  seqnames = "chr3", 
  ranges = IRanges(start = 149910168, end = 149911733)
) # r4
# entire_region.gr <- GRanges(
#   seqnames = "chr13", 
#   ranges = IRanges(start = 42681330, end = 42692763)
# )  # r5
vmrs.gr[countOverlaps(vmrs.gr, entire_region.gr)>0]


