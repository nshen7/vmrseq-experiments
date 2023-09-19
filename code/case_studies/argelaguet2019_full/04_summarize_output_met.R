source("code/SETPATHS.R")
library(BiocParallel)
register(MulticoreParam(workers = 14))

source("code/case_studies/util_functions/4summarizeOutput.R")

read_dir_met <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', 'met')
read_dir_rna <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', 'rna')
write_dir <- here('data', 'interim', 'case_studies', 'argelaguet2019_full', '04_summarize_output_met')
if (!file.exists(write_dir)) dir.create(write_dir)

CHR_NAMES <- as.character(1:19)
upstream <- 2000; downstream <- 2000

# ---- Summarize VMR info ----

# summarizeOutputRegion(read_dir_met, write_dir, methods = c("vmrseq"), chr_names = CHR_NAMES)

# ---- Summarize promoters ----
rna.gr <- granges(readRDS(here(read_dir_rna, 'SingleCellExperiment_rna_vst.rds')))

## Take the [start-upstrem, start-downstream] region as promoters for each gene
rs_promoter.se <- getRegionSummary(
  se_dirs = here(read_dir_met, "vmrseq", "input", CHR_NAMES), 
  regions.gr = promoters(rna.gr, upstream, downstream)
)
saveHDF5SummarizedExperiment(rs_promoter.se, here(write_dir, "promoters_regionSummary"), replace = TRUE)

