source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(argparse)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
n_cores <- 22
register(MulticoreParam(workers = n_cores))

p <- ArgumentParser(description = 'Run vmrseq on Luo 2017 full dataset')
p$add_argument('--chr', type = "integer", help = 'Chromosome number')
p$add_argument('--alpha', type = "double", help = 'Alpha')
args <- p$parse_args(commandArgs(TRUE))
chr <- paste0("chr", args$chr)
alpha <- args$alpha
print(paste0("chr = ", chr, "   alpha = ", alpha))
print(str(alpha))

read_dir <- "data/interim/case_studies/luo2017mice_full/vmrseq/input/"
write_dir <- paste0("data/interim/case_studies/luo2017mice_full_sens_analysis/alpha/alpha", alpha, "/vmrseq/output/")
if (!file.exists(write_dir)) dir.create(write_dir, recursive = T)

# load input
SE <- loadHDF5SummarizedExperiment(paste0(read_dir, chr))

# QC: remove sites with across-cell coverage < 3
total <- rowSums(assays(SE)[[1]]>=0, na.rm = T)
SE <- subset(SE, total >= 3)

# run vmrseq
gr <- vmrseq.smooth(SE)
fit <- vmrseq.fit(gr, alpha = alpha)

# save model output
saveRDS(fit, paste0(write_dir, chr, ".rds"))

