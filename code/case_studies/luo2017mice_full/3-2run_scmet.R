.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/")
source("code/case_studies/util_functions/3-2runScmet.R")

read_dir <- "data/interim/case_studies/luo2017mice_full/scmet/input/scmet_input_Y.txt.gz"
write_dir <- "data/interim/case_studies/luo2017mice_full/scmet/output/"
if (!file.exists(write_dir)) dir.create(write_dir)

plot_dir <- "plots/case_studies/luo2017mice_full/scmet/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr1", n_cores = 22)
runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr2", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr3", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr4", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr5", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr6", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr7", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr8", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr9", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr10", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr11", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr12", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr13", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr14", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr15", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr16", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr17", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr18", n_cores = 34)
# runScmet1Chr(read_dir, write_dir, plot_dir, chr = "chr19", n_cores = 34)