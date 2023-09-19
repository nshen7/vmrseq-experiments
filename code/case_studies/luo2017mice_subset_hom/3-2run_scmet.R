source('code/SETPATHS.R')
source("code/case_studies/util_functions/3-2runScmet.R")

read_dir <- "data/interim/case_studies/luo2017mice_subset_hom/scmet/input/scmet_input_Y.txt.gz"
write_dir <- "data/interim/case_studies/luo2017mice_subset_hom/scmet/output/"
plot_dir <- "plots/case_studies/luo2017mice_subset_hom/scmet/"

runScmet(read_dir, write_dir, plot_dir, efdr = 0.1, n_cores = 28)