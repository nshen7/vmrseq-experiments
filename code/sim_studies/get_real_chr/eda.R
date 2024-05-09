source("code/SETPATHS.R")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(SummarizedExperiment))
register(MulticoreParam(workers = 14))

subtype <- "IT-L23_Cux1"
chromosome <- "chr1"

folder <- paste0("data/interim/sim_studies/real/liu2021_",subtype,"_",chromosome,"/")
pos0 <- fread(paste0(folder, "pos0.txt.gz")) %>% unlist %>% unname
all_file_dirs <- paste0(folder, grep("cell", list.files(folder), value = TRUE))

getCovProp <- function(file_dir) {
  vec <- fread(file_dir) %>% unlist %>% unname
  prop <- sum(!is.na(vec)) / length(vec)
  return(prop)
}

cov_props <- bplapply(all_file_dirs, getCovProp) %>% unlist
df <- data.frame(fileDir = all_file_dirs, cpgCovProp = cov_props)

fwrite(df, paste0("data/interim/sim_studies/real/liu2021_",subtype,"_",chromosome,"_metadata.csv"))


# ==== after running above ====

df <- fread("data/interim/sim_studies/real/liu2021_IT-L23_Cux1_chr1_metadata.csv")

hist(df$cpgCovProp, breaks = 100)

mean(df$cpgCovProp) # 7.4%
sd(df$cpgCovProp) # 2.3%

quantile(df$cpgCovProp, prob = seq(0,1,1/3))
#          0%   33.33333%   66.66667%        100% 
# 0.008881207 0.063019976 0.082245468 0.258942784 

quantile(df$cpgCovProp, prob = seq(0.05,0.95,0.3))
#         5%        35%        65%        95% 
# 0.03807328 0.06403174 0.08125780 0.11252324 

quantile(df$cpgCovProp, prob = seq(0,1,0.1))
#          0%         10%         20%         30%         40%         50% 
# 0.008881207 0.045492400 0.053910476 0.061046336 0.066826228 0.072402599 
#         60%         70%         80%         90%        100% 
# 0.078024009 0.084530795 0.092108956 0.103019399 0.258942784 

sum(df$cpgCovProp >= 0.04 & df$cpgCovProp < 0.06) # 1464
sum(df$cpgCovProp >= 0.06 & df$cpgCovProp < 0.08) # 2253
sum(df$cpgCovProp >= 0.08 & df$cpgCovProp < 0.1) # 1639

sum(df$cpgCovProp >= 0.01 & df$cpgCovProp < 0.05) # 969
sum(df$cpgCovProp >= 0.05 & df$cpgCovProp < 0.09) # 4089
sum(df$cpgCovProp >= 0.09 & df$cpgCovProp < 0.13) # 1406

sum(df$cpgCovProp >= 0.03 & df$cpgCovProp < 0.06) # 1762
sum(df$cpgCovProp >= 0.06 & df$cpgCovProp < 0.09) # 3191
sum(df$cpgCovProp >= 0.09 & df$cpgCovProp < 0.12) # 1310

(qt <- quantile(df$cpgCovProp, prob = seq(0.02,0.98,0.32)))
#         2%        34%        66%        98% 
# 0.03121115 0.06334833 0.08185095 0.12356401 
sum(df$cpgCovProp >= qt[1] & df$cpgCovProp < qt[2]) # 2096
sum(df$cpgCovProp >= qt[2] & df$cpgCovProp < qt[3]) # 2096
sum(df$cpgCovProp >= qt[3] & df$cpgCovProp < qt[4]) # 2096
df %>% filter(cpgCovProp >= qt[1] & cpgCovProp < qt[2]) %>% 
  summarise(mean = mean(cpgCovProp), sd = sd(cpgCovProp)) %>% round(3)
df %>% filter(cpgCovProp >= qt[2] & cpgCovProp < qt[3]) %>% 
  summarise(mean = mean(cpgCovProp), sd = sd(cpgCovProp)) %>% round(3)
df %>% filter(cpgCovProp >= qt[3] & cpgCovProp < qt[4]) %>% 
  summarise(mean = mean(cpgCovProp), sd = sd(cpgCovProp)) %>% round(3)

  