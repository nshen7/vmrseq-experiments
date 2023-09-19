source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(argparse)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
n_cores <- 22
register(MulticoreParam(workers = n_cores))

read_dir <- "data/interim/case_studies/luo2017mice_full/vmrseq/input/"
write_dir <- "data/interim/case_studies/luo2017mice_full/vmrseq/output/"

p <- ArgumentParser(description = 'Formatting for smallwood')
p$add_argument('--chr', type = "integer", help = 'Chromosome number')
args <- p$parse_args(commandArgs(TRUE))
chr <- paste0("chr", args$chr)
print(chr)

# load input
SE <- loadHDF5SummarizedExperiment(paste0(read_dir, chr))

# QC: remove sites with across-cell coverage < 3
total <- rowSums(assays(SE)[[1]]>=0, na.rm = T)
SE <- subset(SE, total >= 3)

# run vmrseq
gr <- vmrseq.smooth(SE)
fit <- vmrseq.fit(gr, alpha = 0.05)

# save model output
saveRDS(fit, paste0(write_dir, chr, ".rds"))

