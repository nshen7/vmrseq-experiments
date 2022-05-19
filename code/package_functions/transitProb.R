# suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(tidyverse))

# ==== class def ====
setClass("transitProbs", representation(max_dist_bp = "numeric",
                                        buffer_bp = "numeric",
                                        transit_probs = "data.frame", 
                                        buffer_probs = "data.frame",
                                        train = "data.frame"))

# ==== helper function ====
.computeProb1Cell <- function(df, max_dist_bp, buffer_bp){
  ### df is a data.frame with 3 columns in strict order of: (chr), (pos), (binary methyl value) 
  
  df <- df %>% na.omit()
  if(!all(df[,3] %in% 0:1)) stop("Methylation value should be either 0 or 1.")

  colnames(df) <- c("chr", "pos", "MF")
  smr_df <- df %>% 
    group_by(chr) %>% 
    mutate(MF_lag1 = lag(MF), dist_bp = c(NA, diff(pos))) %>%
    filter(dist_bp <= max_dist_bp + buffer_bp) %>%
    group_by(dist_bp) %>%
    summarise(N_00 = sum(MF == 0 & MF_lag1==0), # number of sites from 0 to 0
              N_10 = sum(MF == 1 & MF_lag1==0), # number of sites from 0 to 1
              N_01 = sum(MF == 0 & MF_lag1==1), # number of sites from 1 to 0
              N_11 = sum(MF == 1 & MF_lag1==1)  # number of sites from 1 to 1
              ) %>% 
    mutate(p_00 = N_00 / (N_00 + N_10), # P(X_i = 0 | X_{i-1} = 0)
           p_01 = N_01 / (N_01 + N_11) # P(X_i = 0 | X_{i-1} = 1)
           ) %>%
    mutate(p_10 = 1 - p_00, # P(X_i = 1 | X_{i-1} = 0)
           p_11 = 1 - p_01 # P(X_i = 1 | X_{i-1} = 1)
           ) %>%
    select(-c(N_00, N_10, N_01, N_11))
  return(smr_df)
}

.estimTransitProbsFromSummary <- function(smr_cells, max_dist_bp, buffer_bp, cores, ...){
  ### `max_dist_bp` is the maximum CpG-CpG distance that is believed to affect transition probabilities
  ### `buffer_bp` is the length of buffer used to fit smoothing curve near maximum distance
  ### `...` is for additional arguments in loess fitting
  
  stopifnot("max(smr_cells$dist_bp) not equal to max_dist_bp + buffer_bp." = max(smr_cells$dist_bp) == max_dist_bp + buffer_bp)
  
  print("Computing mean and var of the probs across cells.")
  train_data <- smr_cells %>%
    group_by(dist_bp) %>%
    summarise(pbar_00 = mean(p_00, na.rm = T),
              pbar_01 = mean(p_01, na.rm = T),
              pbar_10 = mean(p_10, na.rm = T),
              pbar_11 = mean(p_11, na.rm = T),
              var_00 = var(p_00, na.rm = T), 
              var_01 = var(p_01, na.rm = T),
              var_10 = var(p_10, na.rm = T),
              var_11 = var(p_11, na.rm = T)) %>% 
    add_row(dist_bp = 1, 
            pbar_00 = NA, pbar_01 = NA, pbar_10 = NA, pbar_11 = NA, 
            var_00 = NA, var_01 = NA, var_10 = NA, var_11 = NA,
            .before = 1)
  
  print("Loess smoothing over probs.")
  x <- 1:(max_dist_bp+buffer_bp) # starting from 1 so that row number equal to distance
  smoothed_probs <- with(train_data %>% na.omit() %>% filter(var_00!=0, var_01!=0, var_10!=0, var_11!=0),
                         data.frame(phat_00 = loess(pbar_00 ~ log(dist_bp), weights = 1/var_00, ...) %>% predict(newdata = log(x)),
                                    phat_01 = loess(pbar_01 ~ log(dist_bp), weights = 1/var_01, ...) %>% predict(newdata = log(x)),
                                    phat_10 = loess(pbar_10 ~ log(dist_bp), weights = 1/var_10, ...) %>% predict(newdata = log(x)),
                                    phat_11 = loess(pbar_11 ~ log(dist_bp), weights = 1/var_11, ...) %>% predict(newdata = log(x))
                         )
  )
  tp <- new("transitProbs", 
            max_dist_bp = max_dist_bp, buffer_bp = buffer_bp,
            transit_probs = smoothed_probs[1:max_dist_bp,], 
            buffer_probs = smoothed_probs[(1:buffer_bp)+max_dist_bp,],
            train = train_data %>% select(-dist_bp))
  return(tp)
}

# ==== main ====
### Postfixes rule: P(0|0) => '00'; P(0|1) => '01'; P(1|0) => '10'; P(1|1) => '11'

### Compute the probs conditioning on CpG-CpG distance in each cell and smooth over the probs (each cell has equal weight)
### **rows with NA will be omitted (better input data without NAs to save computing time)
estimTransitProbs <- function(list_cells, max_dist_bp = 2000, buffer_bp = 3000, cores = 1, ...) {
  
  print("Computing transition probs within cells.")
  smr_cells <- do.call(rbind, list_cells %>% mclapply(.computeProb1Cell, 
                                                      mc.cores = cores, 
                                                      max_dist_bp = max_dist_bp, buffer_bp = buffer_bp) 
  ) # proportions of sites in categories 00, 01, 10, 11 conditioning on CpG distance
  
  return(.estimTransitProbsFromSummary(smr_cells = smr_cells, 
                                       max_dist_bp = max_dist_bp, buffer_bp = buffer_bp, 
                                       cores = cores, ...))
}


### Diagnostic plots of transition prob fitting
plotTransitProbs <- function(tp, fitted = T, linewidth = 0.5){
  ### `fitted` is logical value indiating whether to add fitted loess smooth line
  ### `linewidth` is the width of fitted loess smooth line, only applicable when `fitted = T`
  if(nrow(tp@transit_probs) != tp@max_dist_bp) 
    stop("Number of rows in `transit_probs` not equal to `maxdist_bp`.")
  if(nrow(tp@buffer_probs) != tp@buffer_bp) 
    stop("Number of rows in `buffer_probs` not equal to `buffer_bp`.")
  if(nrow(tp@transit_probs) + nrow(tp@buffer_probs) != nrow(tp@train))
    stop("Number of rows in `transit_probs` and `buffer_probs` not adding up to number of rows in `train`.")
  
  plot_df <- data.frame(dist_bp = 1:nrow(tp@train), rbind(tp@transit_probs, tp@buffer_probs), tp@train) %>%
    select(-starts_with("var")) %>% 
    pivot_longer(cols = -1, names_to = c(".value", "type"),
                 names_pattern = "(.*)_(.*)")
  
  type_labs <- c("P(0|0)","P(0|1)","P(1|0)","P(1|1)")
  names(type_labs) <- c("00","01","10","11")
  
  if (fitted) {
    plot_df %>% 
      ggplot() + 
      geom_point(aes(dist_bp, pbar), color = "grey", size = 0.2) + 
      geom_vline(xintercept = tp@max_dist_bp, color = "light blue", linetype = "dashed") +
      geom_path(aes(dist_bp, phat), color = "red", size = 0.2) + 
      theme_bw() + ylim(0, 1) +
      xlab("CpG-CpG Distance") + ylab("Transition Probability") +
      facet_wrap(~ type, labeller = labeller(type = type_labs)) 
  } else {
    plot_df %>% 
      ggplot() + 
      geom_point(aes(dist_bp, pbar), color = "grey", size = 0.2) + 
      geom_vline(xintercept = tp@max_dist_bp, color = "light blue", linetype = "dashed") +
      theme_bw() + ylim(0, 1) +
      xlab("CpG-CpG Distance") + ylab("Transition Probability") +
      facet_wrap(~ type, labeller = labeller(type = type_labs)) 
  }
} 


# # ==== example ====
# cells.se <- loadHDF5SummarizedExperiment(dir = "data/processed/processed_luo2017_mice/subtype_mL23_649cells")
# df1 <- data.frame(chr = seqnames(cells.se),
#                   pos = start(cells.se),
#                   MF = round(assays(cells.se[,1])$M / assays(cells.se[,1])$Cov))
# df2 <- data.frame(chr = seqnames(cells.se),
#                   pos = start(cells.se),
#                   MF = round(assays(cells.se[,2])$M / assays(cells.se[,2])$Cov))
# df3 <- data.frame(chr = seqnames(cells.se),
#                   pos = start(cells.se),
#                   MF = round(assays(cells.se[,3])$M / assays(cells.se[,3])$Cov))
# list_cells <- list(df1, df2, df3)
# 
# tp <- estimTransitProbs(list_cells, cores = 8, span = 0.5, degree = 1)
# 
# ## sanity check
# all.equal(tp@transit_probs[-1,]$phat_00 + tp@transit_probs[-1,]$phat_10, rep(1, nrow(tp@transit_probs[-1,])))
# head(tp@transit_probs)
# 
# write_rds(tp, here::here("code/package_functions/transitProbs_3cells_mL23_Luo2017.rds"))
# 
# plotTransitProbs(tp)
# ggsave(here::here("code/package_functions/transitProbs_3cells_mL23_Luo2017.png"), width = 7, height = 5)
# 
