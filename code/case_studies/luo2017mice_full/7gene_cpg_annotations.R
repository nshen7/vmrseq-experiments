source("code/SETPATHS.R")
library(annotatr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(HDF5Array)

read_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

# marker_genes <- read.csv("data/metadata/metadata_luo2017/luo_supp_tables/NIHMS893063-supplement-Table_S3_csv.csv", skip = 1)
# idx <- grep("Cluster", colnames(marker_genes))
# colnames(marker_genes)[idx] <- marker_genes[1,idx]
# marker_genes <- marker_genes[-1, ]

res.gr <- list(
  'vseq' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_vmrs"))),
  'vseq_cr' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "vmrseq_regionSummary_crs"))),
  'scbs' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "scbs_regionSummary_vmrs"))),
  'smwd' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "smallwood_regionSummary_vmrs"))),
  'scmet' = granges(loadHDF5SummarizedExperiment(paste0(read_dir, "scmet_regionSummary_vmrs")))
)

methodName <- function(method) 
  switch (method,
    'vseq' = 'vmrseq',
    'vseq_cr' = 'vmrseq CRs',
    'scbs' = 'scbs',
    'smwd' = 'Smallwood',
    'scmet' = 'scMET'
  )

# ==== load gene annotations ====

types_gene <- grep("mm10_genes.*", builtin_annotations(), value = T)
annt_genes <- build_annotations(genome = 'mm10', annotations = types_gene)

# ==== annotate detected regions with genes ====

getGeneAnnot <- function(method, top_n_regions){
  
  regions.gr <- res.gr[[method]]
  if (method != "vseq_cr") {
    metric <- switch(method,
                     "vseq" = regions.gr$loglik_diff,
                     "scbs" = regions.gr$mcols.var,
                     "smwd" = regions.gr$var_lb,
                     "scmet" = regions.gr$gamma)
    regions.gr <- regions.gr[order(metric, decreasing = TRUE)]
  }
  
  if (!is.null(top_n_regions)) {
    if (method %in% c("vseq_cr")) stop("Method does not provide rank!")
    regions.gr <- regions.gr[1:top_n_regions]
  }
  
  if (method == "smwd") {
    regions.gr <- GenomicRanges::reduce(regions.gr)
  }
  
  
  method_annt <- annotate_regions(
    regions = regions.gr,
    annotations = annt_genes,
    ignore.strand = TRUE,
    quiet = FALSE
  ) %>% data.frame()
  
  return(list(regions.gr, method_annt))
}

## Plot the distribution of gene annotation types compared with randomized regions
plotGeneAnnotDist <- function(method, top_n_regions = NULL) {
  
  output <- getGeneAnnot(method = method, top_n_regions = top_n_regions)
  regions.gr <- output[[1]]
  method_annt <- output[[2]]
  
  # Generate random regions to check enrichments
  seqinfo(regions.gr) <- Seqinfo(seqnames(seqinfo(regions.gr)), genome = "mm10")
  rand_regions <- randomize_regions(
    regions = regions.gr,
    allow.overlaps = TRUE,
    per.chromosome = TRUE)
  
  # Annotate the random regions using the same annotations as above
  set.seed(100)
  rand_annt <- annotate_regions(
    regions = rand_regions,
    annotations = annt_genes,
    ignore.strand = TRUE,
    quiet = TRUE)
  
  title_seg <- ifelse(
    is.null(top_n_regions), 
    yes = paste0(" - all ", length(regions.gr), " regions"), 
    no = paste0(" - top ", top_n_regions, " regions")
  )
  plot_title <- paste0(methodName(method), title_seg)
  plot_annotation(
    annotated_regions = method_annt,
    annotated_random = rand_annt,
    annotation_order = unique(annt_genes$type),
    plot_title = plot_title,
    x_label = 'Annotations',
    y_label = 'Count')
  seg <- ifelse(is.null(top_n_regions), yes = "", no = paste0("_top", top_n_regions, "regions"))
  ggsave(paste0(plot_dir, "barplot_annot_gene_dist_vs_random_", method, seg, ".png"),
         width = 8, height = 6)
}

plotGeneAnnotDist(method = "vseq")
plotGeneAnnotDist(method = "vseq_cr")
plotGeneAnnotDist(method = "scbs")
plotGeneAnnotDist(method = "smwd")
plotGeneAnnotDist(method = "scmet")

plotGeneAnnotDist(method = "vseq", top_n_regions = 1000)
plotGeneAnnotDist(method = "scbs", top_n_regions = 1000)
plotGeneAnnotDist(method = "smwd", top_n_regions = 1000)
plotGeneAnnotDist(method = "scmet", top_n_regions = 1000)

plotGeneAnnotDist(method = "vseq", top_n_regions = 5000)
plotGeneAnnotDist(method = "scbs", top_n_regions = 5000)
plotGeneAnnotDist(method = "smwd", top_n_regions = 5000)
plotGeneAnnotDist(method = "scmet", top_n_regions = 5000)



# ==== load CpG annotations ====

types_cpg <- "mm10_cpgs"
annt_cpgs <- build_annotations(genome = 'mm10', annotations = types_cpg)

cpg_gene <- annotate_regions(
  regions = annt_cpgs,
  annotations = annt_genes,
  ignore.strand = TRUE,
  quiet = FALSE
) %>% data.frame()


# ==== annotate detected regions with cpgs ====

getCpgAnnot <- function(method, top_n_regions){
  
  regions.gr <- res.gr[[method]]
  if (!method %in% c("vseq_cr")) {
    metric <- switch(method,
                     "vseq" = regions.gr$loglik_diff,
                     "scbs" = regions.gr$mcols.var,
                     "smwd" = regions.gr$var_lb,
                     "scmet" = regions.gr$gamma)
    regions.gr <- regions.gr[order(metric, decreasing = TRUE)]
  }
  
  if (!is.null(top_n_regions)) {
    if (method %in% c("vseq_cr")) stop("Method does not provide rank!")
    regions.gr <- regions.gr[1:top_n_regions]
  }
  
  if (method == "smwd") {
    regions.gr <- GenomicRanges::reduce(regions.gr)
  }
  
  method_annt <- annotate_regions(
    regions = regions.gr,
    annotations = annt_cpgs,
    ignore.strand = TRUE,
    quiet = FALSE
  ) %>% data.frame()
  
  return(list(regions.gr, method_annt))
}

## Plot the distribution of cpg annotation types compared with randomized regions
plotCpgAnnotDist <- function(method, top_n_regions = NULL) {
  
  output <- getCpgAnnot(method = method, top_n_regions = top_n_regions)
  regions.gr <- output[[1]]
  method_annt <- output[[2]]
  
  # Generate random regions to check enrichments
  seqinfo(regions.gr) <- Seqinfo(seqnames(seqinfo(regions.gr)), genome = "mm10")
  rand_regions <- randomize_regions(
    regions = regions.gr,
    allow.overlaps = TRUE,
    per.chromosome = TRUE)
  
  # Annotate the random regions using the same annotations as above
  rand_annt <- annotate_regions(
    regions = rand_regions,
    annotations = annt_cpgs,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  # Plot distribution of annotation types
  title_seg <- ifelse(
    is.null(top_n_regions), 
    yes = paste0(" - all ", length(regions.gr), " regions"), 
    no = paste0(" - top ", top_n_regions, " regions")
  )
  plot_title <- paste0(methodName(method), title_seg)
  plot_annotation(
    annotated_regions = method_annt,
    annotated_random = rand_annt,
    annotation_order = unique(annt_cpgs$type),
    plot_title = plot_title,
    x_label = 'Annotations',
    y_label = 'Count')
  seg <- ifelse(is.null(top_n_regions), yes = "", no = paste0("_top", top_n_regions, "regions"))
  ggsave(paste0(plot_dir, "barplot_annot_cpg_dist_vs_random_", method, seg, ".png"),
         width = 8, height = 6)
}

plotCpgAnnotDist(method = "vseq")
plotCpgAnnotDist(method = "vseq_cr")
plotCpgAnnotDist(method = "scbs")
plotCpgAnnotDist(method = "smwd")
plotCpgAnnotDist(method = "scmet")

plotCpgAnnotDist(method = "vseq", top_n_regions = 1000)
plotCpgAnnotDist(method = "scbs", top_n_regions = 1000)
plotCpgAnnotDist(method = "smwd", top_n_regions = 1000)
plotCpgAnnotDist(method = "scmet", top_n_regions = 1000)

plotCpgAnnotDist(method = "vseq", top_n_regions = 5000)
plotCpgAnnotDist(method = "scbs", top_n_regions = 5000)
plotCpgAnnotDist(method = "smwd", top_n_regions = 5000)
plotCpgAnnotDist(method = "scmet", top_n_regions = 5000)


# ==== annotate detected with both cpgs and genes ====

getBothAnnot <- function(method, top_n_regions){

  regions.gr <- res.gr[[method]]
  if (!method %in% c("vseq_cr")) {
    metric <- switch(method,
                     "vseq" = regions.gr$loglik_diff,
                     "scbs" = regions.gr$mcols.var,
                     "smwd" = regions.gr$var_lb,
                     "scmet" = regions.gr$gamma)
    regions.gr$rank <- order(metric, decreasing = TRUE)
  }

  if (!is.null(top_n_regions)) {
    if (method %in% c("vseq_cr")) stop("Method does not provide rank!")
    top_ind <- regions.gr$rank[1:top_n_regions]
    regions.gr <- regions.gr[top_ind]
  }

  if (method == "smwd") {
    regions.gr <- GenomicRanges::reduce(regions.gr)
  }

  method_annt_cpgs <- annotate_regions(regions = regions.gr,
                                       annotations = annt_cpgs,
                                       ignore.strand = TRUE,
                                       quiet = FALSE) %>%
    data.frame() %>%
    subset_order_tbl(col = "annot.type", col_order = unique(annt_cpgs$type)) %>%
    distinct(across(c("seqnames", "start", "end", "annot.type")),  .keep_all = TRUE)

  method_annt_genes <- annotate_regions(regions = regions.gr,
                                        annotations = annt_genes,
                                        ignore.strand = TRUE,
                                        quiet = FALSE) %>%
    data.frame() %>%
    subset_order_tbl(col = "annot.type", col_order = unique(annt_genes$type)) %>%
    distinct(across(c("seqnames", "start", "end", "annot.type")),  .keep_all = TRUE)

  return(list(regions.gr, method_annt_cpgs, method_annt_genes))
}

## Plot the distribution of cpg annotation types compared with randomized regions
plotBothAnnotDist <- function(method, top_n_regions = NULL) {

  output <- getBothAnnot(method = method, top_n_regions = top_n_regions)
  regions.gr        <- output[[1]]
  method_annt_cpgs  <- output[[2]]
  method_annt_genes <- output[[3]]
  method_annt <- method_annt_cpgs %>%
    left_join(method_annt_genes,
              by = c("seqnames", "start", "end"),
              multiple = "all",
              suffix = c(".cpgs", ".genes")) %>%
    distinct(across(c("seqnames", "start", "end", "annot.type.cpgs", "annot.type.genes")), .keep_all = TRUE)


  # Generate random regions to check enrichments
  seqinfo(regions.gr) <- Seqinfo(seqnames(seqinfo(regions.gr)), genome = "mm10")
  rand_regions <- randomize_regions(
    regions = regions.gr,
    allow.overlaps = TRUE,
    per.chromosome = TRUE)
  # Annotate the random regions using the same annotations as above
  rand_annt_cpgs <- annotate_regions(
    regions = rand_regions,
    annotations = annt_cpgs,
    ignore.strand = TRUE,
    quiet = FALSE) %>%
    as.data.frame() %>%
    subset_order_tbl(col = "annot.type", col_order = unique(annt_cpgs$type))
  rand_annt_genes <- annotate_regions(
    regions = rand_regions,
    annotations = annt_genes,
    ignore.strand = TRUE,
    quiet = FALSE) %>%
    as.data.frame() %>%
    subset_order_tbl(col = "annot.type", col_order = unique(annt_genes$type))
  rand_annt <- rand_annt_cpgs %>%
    left_join(rand_annt_genes,
              by = c("seqnames", "start", "end"),
              multiple = "all",
              suffix = c(".cpgs", ".genes")) %>%
    distinct(across(c("seqnames", "start", "end", "annot.type.cpgs", "annot.type.genes")),  .keep_all = TRUE)

  all_annt <- bind_rows(Data = method_annt, `Random Regions` = rand_annt, .id = "data_type")

  # all_annt %>%
  #   ggplot(aes(x = annot.type.genes, fill = data_type)) +
  #   geom_bar(position = "dodge") +
  #   # geom_histogram(aes(y = after_stat(density))) +
  #   facet_wrap(~ annot.type.cpgs) +
  #   scale_fill_grey() +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 30,hjust = 1),
  #         legend.title = element_blank(),
  #         legend.position = "bottom",
  #         legend.key = element_rect(color = "white"))
  
  title_seg <- ifelse(
    is.null(top_n_regions), 
    yes = paste0(" - all ", length(regions.gr), " regions"), 
    no = paste0(" - top ", top_n_regions, " regions")
  )
  plot_title <- paste0(methodName(method), title_seg)
  props <- all_annt %>%
    group_by(data_type, annot.type.genes, annot.type.cpgs) %>%
    summarise(totals = n()) %>%
    ungroup(annot.type.cpgs) %>%
    summarise(annot.type.cpgs = annot.type.cpgs,
              props.type.cpgs = totals/sum(totals))
  props %>%
    ggplot(aes(x = annot.type.genes, fill = annot.type.cpgs)) +
    geom_col(aes(y = props.type.cpgs), position = "stack") +
    # geom_histogram(aes(y = after_stat(density))) +
    facet_wrap(~ data_type, nrow = 2) +
    # scale_fill_grey() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30,hjust = 1),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.key = element_rect(color = "white")) +
    ylab("Proportion of CpG annotation types") +
    xlab("Gene annotation types") +
    ggtitle(plot_title)
  seg <- ifelse(is.null(top_n_regions), yes = "", no = paste0("_top", top_n_regions, "regions"))
  ggsave(paste0(plot_dir, "barplot_annot_cpg&gene_dist_vs_random_", method, seg, ".png"),
         width = 8, height = 8)
}


# plotBothAnnotDist(method = "vseq")
# plotBothAnnotDist(method = "vseq_cr")
# plotBothAnnotDist(method = "scbs")
# plotBothAnnotDist(method = "smwd")
# plotBothAnnotDist(method = "scmet")
# 
# plotBothAnnotDist(method = "vseq", top_n_regions = 1000)
# plotBothAnnotDist(method = "scbs", top_n_regions = 1000)
# plotBothAnnotDist(method = "smwd", top_n_regions = 1000)
# plotBothAnnotDist(method = "scmet", top_n_regions = 1000)
# 
# plotBothAnnotDist(method = "vseq", top_n_regions = 5000)
# plotBothAnnotDist(method = "scbs", top_n_regions = 5000)
# plotBothAnnotDist(method = "smwd", top_n_regions = 5000)
# plotBothAnnotDist(method = "scmet", top_n_regions = 5000)

 