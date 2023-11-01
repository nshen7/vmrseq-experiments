source("code/SETPATHS.R")
library(BiocParallel)
register(MulticoreParam(workers = 14))

source("code/case_studies/util_functions/4summarizeOutput.R")

read_dir_met <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', 'met')
read_dir_rna <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', 'rna')
write_dir <- here('data', 'interim', 'case_studies', 'argelaguet2019_subset', '04_summarize_output_met')
if (!file.exists(write_dir)) dir.create(write_dir)

CHR_NAMES <- as.character(1:19)

# ---- Summarize VMR info ----

# summarizeOutputRegion(read_dir_met, write_dir, methods = c("vmrseq"), chr_names = CHR_NAMES)


# ---- Summarize promoters ----
library(annotatr)
rna.gr <- granges(readRDS(here(read_dir_rna, 'SingleCellExperiment_rna.rds')))

## This part need to be run in Singularity since network is needed!!!
#### START
# mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
# promoter.gr <- getBiomartCanonicalPromoters(mart = mart, upstream = 2000, downstream = 2000)
# seqlevelsStyle(promoter.gr) <- 'NCBI'
# saveRDS(promoter.gr, here(write_dir, 'GRanges_promoter_info.rds'))
###### END

promoter.gr <- readRDS(here(write_dir, 'GRanges_promoter_info.rds'))
rs_promoter.se <- getRegionSummary(
  se_dirs = here(read_dir_met, "vmrseq", "input", CHR_NAMES), 
  regions.gr = promoter.gr
)
saveHDF5SummarizedExperiment(rs_promoter.se, here(write_dir, "promoters_regionSummary"), replace = TRUE)

