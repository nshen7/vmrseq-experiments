source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(HDF5Array)
library(DelayedMatrixStats)
library(SummarizedExperiment)
library(scales)

read_dir_met <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', '04_summarize_output_met')
read_dir_rna <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', 'rna')
write_dir <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', '06_met_rna_corr')
if (!file.exists(write_dir)) dir.create(write_dir)
plot_dir <- here('plots', 'case_studies', 'argelaguet2019_subset', '06_met_rna_corr')
if (!file.exists(plot_dir)) dir.create(plot_dir)

# ---- Load data ----

md <- fread(here('data', 'metadata', 'argelaguet2019', 'argelaguet2019_subset_met&rna_sample_metadata_processed.csv'))

vmr.se <- loadHDF5SummarizedExperiment(here(read_dir_met, "vmrseq_regionSummary_vmrs"))
colnames(vmr.se) <- md$sample

promoter.se <- loadHDF5SummarizedExperiment(here(read_dir_met, "promoters_regionSummary"))
colnames(promoter.se) <- md$sample
dim(promoter.se) # 20515   726

rna.se <- readRDS(here(read_dir_rna, 'SingleCellExperiment_rna.rds'))
dim(rna.se) # 17265   726

## Subset to genes with present promoter methylation
rna.se <- rna.se[match(granges(promoter.se)$external_gene_name, granges(rna.se)$symbol) %>% na.omit, ] 
dim(rna.se) # 16026   726
all(granges(rna.se)$symbol == granges(promoter.se)$symbol) # TRUE

## Make colnames of vmr.se, promoter.se and rna.se align
all(colnames(rna.se) == md$id_rna) # TRUE
colnames(rna.se) <- md$sample 


quantile(width(granges(vmr.se)))
# 0%   25%   50%   75%  100% 
# 15   732  1166  1768 22126 
quantile(width(granges(promoter.se)))
#   0%  25%  50%  75% 100% 
# 4000 4000 4000 4000 4000 

sum(countOverlaps(granges(vmr.se), granges(promoter.se)) > 0) / nrow(vmr.se) 
# ~ 4.2% VMRs overlap with promoters
sum(countOverlaps(granges(promoter.se), granges(vmr.se)) > 0) / nrow(promoter.se) 
# ~ 6.7% promoters overlap with VMRs

## Mark the VMRs that overlap with promoters
rowData(vmr.se)$olap_promoter <- countOverlaps(granges(vmr.se), granges(promoter.se)) > 0

# ---- VMRs as ref: Compare VMRs and promoters in terms of correlation with gene expression ----

## Use VMRs as reference: pick the closest gene to each VMR
near_vmr.hits <- distanceToNearest(granges(vmr.se), granges(rna.se)) %>%
  as.data.frame() %>%
  filter(distance <= max_dist)
vmr_near.se      <- vmr.se     [near_vmr.hits$queryHits, ]
rna_near.se      <- rna.se     [near_vmr.hits$subjectHits, ]
promoter_near.se <- promoter.se[near_vmr.hits$subjectHits, ]

vmr.mtx <- (assays(vmr_near.se)$M / assays(vmr_near.se)$Cov) %>% as.matrix
rna.mtx <- (assays(rna_near.se)$logcounts) %>% as.matrix
promoter.mtx <- (assays(promoter_near.se)$M / assays(promoter_near.se)$Cov) %>% as.matrix

N <- nrow(vmr.mtx)
corr_vmr      <- map_dbl(.x = 1:N, 
                         .f = ~ cor(vmr.mtx[.x, ], rna.mtx[.x, ], method = 'spearman', use = 'pairwise.complete.obs'))
corr_promoter <- map_dbl(.x = 1:N, 
                         .f = ~ cor(promoter.mtx[.x, ], rna.mtx[.x, ], method = corr_method, use = 'pairwise.complete.obs'))

n_avail_cell_vmr      <- map_dbl(.x = 1:N, 
                                 .f = ~ sum(!is.na(vmr.mtx[.x, ])      & !is.na(rna.mtx[.x, ])))
n_avail_cell_promoter <- map_dbl(.x = 1:N, 
                                 .f = ~ sum(!is.na(promoter.mtx[.x, ]) & !is.na(rna.mtx[.x, ])))

rowData(vmr_near.se) <- cbind(
  rowData(vmr_near.se),
  data.frame(
    rna_gene_name          = granges(rna_near.se)$symbol,
    rna_gene_start         = start(rna_near.se),
    rna_gene_end           = end(rna_near.se),
    rna_gene_strand        = strand(rna_near.se),
    vmr_to_gene_dist       = near_vmr.hits$distance,
    vmr_rna_corr           = corr_vmr, 
    vmr_n_avail_cell       = n_avail_cell_vmr,
    loglik_diff            = granges(vmr_near.se)$loglik_diff,
    promoter_rna_corr      = corr_promoter,
    promoter_n_avail_cell  = n_avail_cell_promoter
  )
)

saveRDS(granges(vmr_near.se), here(write_dir, paste0('GRanges_correlation_', corr_method, '_metNrna_vmrsAsRef_maxDist', max_dist, 'bp.rds')))

for (top_n_regions in c(1000, 'all')) {
  
  vmr_corr.df <- readRDS(here(write_dir, paste0('GRanges_correlation_', corr_method, '_metNrna_vmrsAsRef_maxDist', max_dist, 'bp.rds'))) %>% 
    as.data.frame %>%
    arrange(desc(loglik_diff))
  
  name_seg      <- ifelse(top_n_regions == 'all', yes = '', no = paste0('_top', top_n_regions, 'regions'))
  top_n_regions <- ifelse(top_n_regions == 'all', yes = nrow(vmr_corr.df), no = as.integer(top_n_regions))
  
  vmr_corr.df[1:top_n_regions, ] %>%
    filter(vmr_n_avail_cell >= 10 & promoter_n_avail_cell >= 10) %>% ## ensure min 10 available cells covered
    ggplot(aes(vmr_rna_corr, promoter_rna_corr)) +
    geom_rug(col = 'darkseagreen', alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') +
    geom_vline(xintercept = 0, color = 'darkgrey', linetype = 'dashed') +
    geom_hline(yintercept = 0, color = 'darkgrey', linetype = 'dashed') +
    geom_point(size = 0.8, alpha = 0.5, color = 'cornsilk4') +
    # geom_point(aes(color = olap_promoter), size = 0.8, alpha = 0.5) +
    geom_smooth(method = 'lm', se = F, color = 'darkolivegreen') +
    xlim(-1, 1) + ylim(-1, 1) +
    ggtitle(paste0('# VMRs: ', top_n_regions, '; Max dist.: ', max_dist, ' bp')) + 
    xlab('Correlation of gene expression with nearest VMR') +
    ylab('Correlation of gene expression with promoter') +
    theme_classic() 
  ggsave(here(plot_dir, paste0('correlation_', corr_method, '_metNrna_vmrsAsRef', name_seg, '_maxDist', max_dist, 'bp.png')), width = 5, height = 5)
}

# ---- genes as ref: Compare VMRs and promoters in terms of correlation with gene expression ----

## Use genes as reference: pick the closest VMR to each gene
## (According to ?nearest, are there multiple overlaps, one of the overlapping ranges is chosen arbitrarily.)
near_rna.hits <- distanceToNearest(granges(rna.se), granges(vmr.se)) %>%
  as.data.frame() %>%
  filter(distance <= max_dist)
vmr_near.se      <- vmr.se     [near_rna.hits$subjectHits, ]
promoter_near.se <- promoter.se[near_rna.hits$queryHits, ]
rna_near.se      <- rna.se     [near_rna.hits$queryHits, ]

vmr.mtx      <- (assays(vmr_near.se)$M      / assays(vmr_near.se)$Cov) %>% as.matrix
promoter.mtx <- (assays(promoter_near.se)$M / assays(promoter_near.se)$Cov) %>% as.matrix
rna.mtx      <- (assays(rna_near.se)$logcounts) %>% as.matrix

## Compute two types of correlation between RNA and methylation, RNA and promoter
N <- nrow(rna.mtx)
corr_vmr      <- map_dbl(.x = 1:N, 
                         .f = ~ cor(vmr.mtx[.x, ],      rna.mtx[.x, ], method = corr_method, use = 'pairwise.complete.obs'))
corr_promoter <- map_dbl(.x = 1:N, 
                         .f = ~ cor(promoter.mtx[.x, ], rna.mtx[.x, ], method = corr_method, use = 'pairwise.complete.obs'))

n_avail_cell_vmr      <- map_dbl(.x = 1:N, 
                                 .f = ~ sum(!is.na(vmr.mtx[.x, ])      & !is.na(rna.mtx[.x, ])))
n_avail_cell_promoter <- map_dbl(.x = 1:N, 
                                 .f = ~ sum(!is.na(promoter.mtx[.x, ]) & !is.na(rna.mtx[.x, ])))

rowData(vmr_near.se) <- cbind(
  rowData(vmr_near.se),
  data.frame(
    rna_gene_name          = granges(rna_near.se)$symbol,
    rna_gene_start         = start(rna_near.se),
    rna_gene_end           = end(rna_near.se),
    rna_gene_strand        = strand(rna_near.se),
    vmr_to_gene_dist       = near_rna.hits$distance,
    vmr_rna_corr           = corr_vmr, 
    vmr_n_avail_cell       = n_avail_cell_vmr,
    loglik_diff            = granges(vmr_near.se)$loglik_diff,
    promoter_rna_corr      = corr_promoter,
    promoter_n_avail_cell  = n_avail_cell_promoter
  )
)
saveRDS(granges(vmr_near.se), here(write_dir, paste0('GRanges_correlation_', corr_method, '_metNrna_genesAsRef_maxDist', max_dist, 'bp.rds')))

## Plotting
for (top_n_regions in c(1000, 'all')) {
  
  vmr_corr.df <- readRDS(here(write_dir, paste0('GRanges_correlation_', corr_method, '_metNrna_genesAsRef_maxDist', max_dist, 'bp.rds'))) %>% 
    as.data.frame %>%
    arrange(desc(loglik_diff))
  
  name_seg      <- ifelse(top_n_regions == 'all', yes = '', no = paste0('_top', top_n_regions, 'regions'))
  top_n_regions <- ifelse(top_n_regions == 'all', yes = nrow(vmr_corr.df), no = as.integer(top_n_regions))
  
  vmr_corr.df[1:top_n_regions, ] %>%
    filter(vmr_n_avail_cell >= 10 & promoter_n_avail_cell >= 10) %>% ## ensure min 10 available cells covered
    ggplot(aes(vmr_rna_corr, promoter_rna_corr)) +
    geom_rug(col = 'darkseagreen', alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') +
    geom_vline(xintercept = 0, color = 'darkgrey', linetype = 'dashed') +
    geom_hline(yintercept = 0, color = 'darkgrey', linetype = 'dashed') +
    geom_point(size = 0.8, alpha = 0.5, color = 'cornsilk4') +
    geom_smooth(method = 'lm', se = F, color = 'darkolivegreen') +
    xlim(-1, 1) + ylim(-1, 1) +
    ggtitle(paste0('# genes: ', top_n_regions, '; Max dist.: ', max_dist, ' bp')) + 
    xlab('Correlation of gene expression with nearest VMR') +
    ylab('Correlation of gene expression with promoter') +
    theme_classic() 
  ggsave(here(plot_dir, paste0('correlation_', corr_method, '_metNrna_genesAsRef', name_seg, '_maxDist', max_dist, 'bp.png')), width = 5, height = 5)
  
}

