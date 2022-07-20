.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
library(GenomicRanges)
library(bumphunter)
library(tidyverse)
library(BiocParallel)
library(vmrseq)

N <- 1000
chr <- sort(sample(1:5, N, replace=TRUE))
pos <- cumsum(round(rexp(N, rate = 1/800) %>% pmax(2)))
# pos <- cumsum(round(rexp(N, rate = 1/800)))

gr <- GRanges(seqnames = Rle(factor(chr)), ranges = IRanges(start = pos, end = pos))
gr$meth <- sample(8:10, length(gr), replace = T)
gr$total <- sample(10:20, length(gr), replace = T)

head(gr)

system.time(vmrseq(gr, BPPARAM = MulticoreParam(workers = 4)))
system.time(vmrseq(gr, BPPARAM = SnowParam(workers = 1)))

# gr %>% as.data.frame() %>% head %>% write.table("gr_head.txt")
#
# minCov = 5;
# cutoff = 0.1; # param for CR calling
# maxGap = 1000; minNumRegion = 5; # params for VMR calling
# smooth = TRUE; maxGapSmooth = 2500; # params for smoother
# bpSpan = 1000; minInSpan = 10; # params for smoother
# tp = NULL;
# maxNumMerge = 1; minNumLong = 10;
# control = vmrseq.control();
# verbose = TRUE
# # BPPARAM = bpparam()
# BPPARAM = SnowParam(workers = 6, type = "FORK")
# # BPPARAM = SnowParam(workers = 1, type = "FORK")
#
# system.time(bplapply(1:10, sqrt))

# system.time(bplapply(1:100000, sqrt, BPPARAM = DoparParam()))
# system.time(bplapply(1:100000, sqrt, BPPARAM = MulticoreParam(workers = 6)))
# system.time(lapply(1:100000, sqrt))
# system.time(parallel::mclapply(1:100000, sqrt, mc.cores = 6))
