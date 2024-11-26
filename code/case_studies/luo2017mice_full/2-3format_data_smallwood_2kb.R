source("code/SETPATHS.R")
library(argparse)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
register(MulticoreParam(workers = 22))

p <- ArgumentParser(description = 'Formatting for smallwood')
p$add_argument('--chr', type = "integer", help = 'Chromosome number')
args <- p$parse_args(commandArgs(TRUE))
chr <- paste0("chr", args$chr)
print(chr)

bp_window <- 2000
bp_step <- 600

write_dir <- "data/interim/case_studies/luo2017mice_full/smallwood_2kb/input/"
if (!file.exists(write_dir)) dir.create(write_dir, recursive = TRUE)
read_dir <- "data/interim/case_studies/luo2017mice_full/vmrseq/input/"
# Load raw data
se <- loadHDF5SummarizedExperiment(dir = paste0(read_dir, chr))

# Format into smallwood input
cuts <- seq(start(se)[1], start(se)[length(se)], bp_step)
len <- length(cuts)
cuts <- cuts[-((len-3):len)]
wds.gr <- GRanges(
  seqnames = seqnames(se)[1],
  ranges = IRanges(start = cuts,
                   end = c(cuts[-length(cuts)]+bp_window, start(se)[length(se)]))
)
hits <- findOverlaps(granges(se), wds.gr)
Indexes <- lapply(
  unique(subjectHits(hits)), 
  function(i) queryHits(hits)[subjectHits(hits)==i]
)
print("Finished finding hits.")
print(Sys.time())

# Divide M_mat into groups to read into RAM and save time
grp_size <- 200000
n_grp <- length(Indexes) %/% grp_size
if (length(Indexes) %% grp_size > 0) n_grp <- n_grp + 1

M <- NULL; Cov <- NULL
for (i in 1:n_grp) {
  grp_idx <- (grp_size*(i-1)+1):min(grp_size*i, length(Indexes))
  head <- Indexes[[grp_idx[1]]]
  tail <- Indexes[[grp_idx[length(grp_idx)]]]
  site_idx <- head[1]:tail[length(tail)]
  grp.se <- se[site_idx,]
  grp_mat <- assays(grp.se)$M_mat %>% as("sparseMatrix")
  
  getFeature <- function(j, type) { # i th feature/window
    mat <- matrix(grp_mat[Indexes[[j]]-site_idx[1]+1, ], ncol = ncol(se))
    if (type == "M") return(round(colSums(mat))) 
    else if (type == "Cov") return(colSums(mat > 0))
    else stop("Wrong 'type' value. Either 'Cov' or 'M'.")
  }
  M <- rbind(M, do.call(
    rbind,
    bplapply(grp_idx, getFeature, type = "M")
  ))
  Cov <- rbind(Cov, do.call(
    rbind,
    bplapply(grp_idx, getFeature, type = "Cov")
  ))
}
print("Finished computing features.")
print(Sys.time())

saveHDF5SummarizedExperiment(
  x = SummarizedExperiment(
    assays = list("M" = M, "Cov" = Cov),
    rowRanges = wds.gr[unique(subjectHits(hits))]
  ),
  dir = paste0(write_dir, chr),
  replace = TRUE
)


