source("code/SETPATHS.R")
library(annotatr)
library(HDF5Array)
library(biomaRt)
library(SummarizedExperiment)
library(GenomicRanges)

read_dir_1 <- "data/metadata/histone_peaks/"
read_dir_2 <- "data/interim/case_studies/luo2017mice_full/result_summary/"
write_dir <- "data/interim/case_studies/luo2017mice_full/result_summary/"
plot_dir <- "plots/case_studies/luo2017mice_full/comparison/"
if (!file.exists(plot_dir)) dir.create(plot_dir)

upstream <- 2000; downstream <- 2000
dissim_metric <- 'manhattan'
res.se <- loadHDF5SummarizedExperiment(here(read_dir_2, "vmrseq_regionSummary_vmrs"))


# ---- RUN THIS WITH INTERNET ACCESS: get TSS info ----
# > listEnsemblArchives()
#              name     date                                 url version
# 1  Ensembl GRCh37 Feb 2014          https://grch37.ensembl.org  GRCh37
# 2     Ensembl 112 May 2024 https://may2024.archive.ensembl.org     112
# 3     Ensembl 111 Jan 2024 https://jan2024.archive.ensembl.org     111
# 4     Ensembl 110 Jul 2023 https://jul2023.archive.ensembl.org     110
# 5     Ensembl 109 Feb 2023 https://feb2023.archive.ensembl.org     109
# 6     Ensembl 108 Oct 2022 https://oct2022.archive.ensembl.org     108
# 7     Ensembl 107 Jul 2022 https://jul2022.archive.ensembl.org     107
# 8     Ensembl 106 Apr 2022 https://apr2022.archive.ensembl.org     106
# 9     Ensembl 105 Dec 2021 https://dec2021.archive.ensembl.org     105
# 10    Ensembl 104 May 2021 https://may2021.archive.ensembl.org     104
# 11    Ensembl 103 Feb 2021 https://feb2021.archive.ensembl.org     103
# 12    Ensembl 102 Nov 2020 https://nov2020.archive.ensembl.org     102
# 13    Ensembl 101 Aug 2020 https://aug2020.archive.ensembl.org     101
# 14    Ensembl 100 Apr 2020 https://apr2020.archive.ensembl.org     100
# 15     Ensembl 99 Jan 2020 https://jan2020.archive.ensembl.org      99
# 16     Ensembl 98 Sep 2019 https://sep2019.archive.ensembl.org      98
# 17     Ensembl 97 Jul 2019 https://jul2019.archive.ensembl.org      97
# 18     Ensembl 80 May 2015 https://may2015.archive.ensembl.org      80
# 19     Ensembl 77 Oct 2014 https://oct2014.archive.ensembl.org      77
# 20     Ensembl 75 Feb 2014 https://feb2014.archive.ensembl.org      75
# 21     Ensembl 54 May 2009 https://may2009.archive.ensembl.org      54

# ensembl80 <- useMart("ENSEMBL_MART_ENSEMBL", 
#                      host = "may2015.archive.ensembl.org",
#                      dataset = "mmusculus_gene_ensembl")
# 
# listAttributes(ensembl80)
# # listFilters(ensembl80)
# 
# attributes <- c("ensembl_gene_id", "external_gene_name", "chromosome_name", "transcription_start_site")
# tss_info <- getBM(attributes = attributes,
#                    mart = ensembl80)
# fwrite(tss_info, here(read_dir_1, "TSS_mmusculus_gene_ensembl_version80.csv"))



# ---- compute NN score for VMRs overlapping with targeted regions ----


getTargetedRegions <- function(target) {
  if (target == 'promoter') {
    tss_info.gr <- fread(here(read_dir_1, "TSS_mmusculus_gene_ensembl_version80.csv")) |> 
      as.data.frame() |> 
      makeGRangesFromDataFrame(start.field="transcription_start_site", end.field = "transcription_start_site")
    seqlevelsStyle(tss_info.gr) <- "UCSC"
    target.gr <- promoters(tss_info.gr, upstream, downstream)
  } else if (target == 'H3K27ac') {
    target.gr <- rtracklayer::import(here(read_dir_1, "ENCFF314HSG.bigBed"))
  } else if (target == 'H3K4me1') {
    target.gr <- rtracklayer::import(here(read_dir_1, "ENCFF975QSJ.bigBed"))
  }
  
  return(target.gr)
}

computeDissimMatrix <- function(target) {
  
  target.gr <- getTargetedRegions(target)
  se <- res.se[countOverlaps(granges(res.se), target.gr) > 0, ]
  MF <- as.matrix(assays(se)$M/assays(se)$Cov) 
  
  # Get dissimilarity matrix for columns (i.e., cells)
  d_mat <- cluster::daisy(t(MF), metric = dissim_metric, stand = FALSE) %>% as.matrix()
  
  fwrite(d_mat, here(write_dir, paste0("dissimilarity_matrix_regional_methyl_vmrseq_", target, ".txt.gz")),
         col.names = F, row.names = F, quote = F)
}

# computeDissimMatrix('promoter')
computeDissimMatrix('H3K27ac')
computeDissimMatrix('H3K4me1')

# ---- compute NN score on VMRs overlapping with targeted reginos ----

source("code/case_studies/util_functions/6nnScore.R")
md <- fread("data/metadata/metadata_luo2017/sample_info_processed.csv")
k <- 100; theta <- 0.7

computeNNScore <- function(target) {
  
  # Count VMRs
  target.gr <- getTargetedRegions(target)
  NVMRs <- length(findOverlaps(granges(res.se), target.gr))
  
  # Compute NN score
  d_mat <- fread(here(read_dir_2, paste0("dissimilarity_matrix_regional_methyl_vmrseq_", target, ".txt.gz")))
  NNscore_broad <- nnScore(d_mat, true_clust = md$Neuron_type1, k, theta)
  NNscore_sub   <- nnScore(d_mat, true_clust = md$Neuron.type,  k, theta)
  
  return(c(NVMRs, NNscore_broad, NNscore_sub))
} 

score.df <- data.frame(
  target = c('promoter', 'H3K27ac', 'H3K4me1'),
  NVMRs = 0,
  NNScore_broad = 0,
  NNScore_sub = 0
)

for (i in 1:nrow(score.df)) {
  score.df[i, 2:4] <- computeNNScore(score.df$target[i])
  print(i)
}
fwrite(score.df, here(write_dir, paste0("nearest_neighbor_score_k", k, "_theta", theta, "_promoter_and_histone.csv")))

# ---- results ----

score_full_broad <- fread(here(write_dir, paste0("nearest_neighbor_score_broadCellType_k", k, "_theta", theta, ".csv"))) %>%
  filter(Method == 'vmrseq') %>%
  filter(NTopRegions == max(NTopRegions))
score_full_sub   <- fread(here(write_dir, paste0("nearest_neighbor_score_subCellType_k", k, "_theta", theta, ".csv"))) %>%
  filter(Method == 'vmrseq') %>%
  filter(NTopRegions == max(NTopRegions)) 

score.df <- fread(here(write_dir, paste0("nearest_neighbor_score_k", k, "_theta", theta, "_promoter_and_histone.csv")))
score.df <- score.df %>% rbind(data.frame(target = 'all VMRs',
                                          NVMRs = score_full_broad$NTopRegions,
                                          NNScore_broad = score_full_broad$NNScore,
                                          NNScore_sub = score_full_sub$NNScore))

score.df[c(4,3,2,1),]
#      target NVMRs NNScore_broad NNScore_sub
# 1: all VMRs 34017     0.9973933   0.7184751
# 2:  H3K4me1 10688     0.9990225   0.6757902
# 3:  H3K27ac  4961     0.9996742   0.6226784
# 4: promoter  6046     0.9980450   0.5457804


