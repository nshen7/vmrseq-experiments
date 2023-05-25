source("code/SETPATHS.R")
library(argparse)
source("code/case_studies/util_functions/3-2runScmet.R")

read_dir <- "data/interim/case_studies/luo2017mice_full/scmet/input/scmet_input_Y.txt.gz"
write_dir <- "data/interim/case_studies/luo2017mice_full/scmet/output/"
if (!file.exists(write_dir)) dir.create(write_dir)

plot_dir <- "plots/case_studies/luo2017mice_full/scmet/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

p <- ArgumentParser(description = 'Run scMET on Luo 2017 full dataset')
p$add_argument('--chr', type = "integer", help = 'Chromosome number')
args <- p$parse_args(commandArgs(TRUE))
chr <- paste0("chr", args$chr)

runScmet1Chr(read_dir = read_dir, write_dir = write_dir, plot_dir, chr = chr, n_cores = 22)

# Combine the HVFs found from each chr
hvfScmet(read_dir = write_dir, write_dir = write_dir)
