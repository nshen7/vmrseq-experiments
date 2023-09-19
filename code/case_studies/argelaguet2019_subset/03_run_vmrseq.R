source("code/SETPATHS.R")
devtools::load_all("../vmrseq-package/vmrseq/")
library(argparse)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
n_cores <- 22
register(MulticoreParam(workers = n_cores))

# ---- Bash arguments ----
p <- ArgumentParser(description = 'Run vmrseq on Argelaguet2019 subset dataset')
p$add_argument('--chr', type = "integer", help = 'Chromosome number')
args <- p$parse_args(commandArgs(TRUE))
chr <- as.character(args$chr)

# ---- main ----
read_dir <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', 'met', 'vmrseq', 'input')
write_dir <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', 'met', 'vmrseq', 'output')
if (!file.exists(write_dir)) dir.create(write_dir)

# load input
SE <- loadHDF5SummarizedExperiment(here(read_dir, chr))
print(dim(SE))

# QC: remove sites with across-cell coverage < 3
total <- rowSums(assays(SE)[[1]] > 0, na.rm = T)
print(table(total))
SE <- subset(SE, total >= 3)

# run vmrseq
gr <- vmrseq.smooth(SE)
fit <- vmrseq.fit(gr, alpha = 0.05)

# save model output
saveRDS(fit, here(write_dir, paste0(chr, ".rds")))
