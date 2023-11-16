.libPaths("/home/nshen7/R/rstudio_4_2_0-biocon_3_15")
setwd("/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/") # path folder is vmrseq-experiments
source("code/sim_studies/helper_functions/(deprecated)helper_functions_generate_chr.R")
devtools::load_all("../vmrseq-package/vmrseq/")

my_eval <- function(res_gr) {
  cat("Site-level power =", with(res_gr, sum(is_vml&!is.na(vmr_index))/sum(is_vml)), "; ")
  cat("Site-level FDR =", with(res_gr, sum(!is_vml&!is.na(vmr_index))/sum(!is.na(vmr_index))), "\n")
  
  smr_true <- res_gr %>% data.frame %>% 
    filter(is_vml) %>%
    group_by(vmr_name) %>%
    summarise(n_cpg = n(), 
              n_dectected = sum(!is.na(vmr_index)), 
              n_incr = sum(!is.na(cr_index)),
              cr_name = max(cr_index, na.rm = T),
              pi1 = max(pi1, na.rm = T),
              med_total = median(total),
              loglik_diff = max(loglik_diff, na.rm = T))
  cat("Region-level power =", sum(smr_true$n_dectected > 0) / nrow(smr_true), "; ")
  smr_detected <- res_gr %>% data.frame %>% 
    filter(!is.na(vmr_index)) %>%
    group_by(vmr_index) %>%
    summarise(cr_index = max(cr_index),
              n_cpg = n(), 
              n_true = sum(!is.na(vmr_name)), 
              pi1 = max(pi1, na.rm = T),
              optim_pi = max(optim_pi),
              mean_MF = mean(meth/total),
              loglik_diff = max(loglik_diff))
  cat("Region-level FDR =", 1 - sum(smr_detected$n_true > 0) / nrow(smr_detected), "\n")
  return(list(true = smr_true, detected = smr_detected))
}

# N <- 500
# cat("N =", N, "\n")
# NP <- 10
# cat("NP =", NP, "\n")
# NV = 100
# cat("NV =", NV, "\n\n")

# for (N in c(100,500,1000,2000)) {
for (N in c(2000)) {
  cat("N =", N, "\n")
# for (r in c(1, 5, 10, 20)) {
  # cat("r =", r, "\n")
  
  # se <- generatePseudoChr(N = N, NV = NV, NP = NP, seed = 2020, nb_r_1g = r)
  se <- generateNullPseudoChr(N = N, NV = 1000, seed = 2022)
  gr <- granges(se)
  gr$meth <- rowSums(assays(se)$MF, na.rm = T)
  gr$total <- rowSums(assays(se)$MF >= 0, na.rm = T)
  
  res_vmrseq <- vmrseq(gr, penalty = 0, tp = readRDS("code/estim_transitProb/tp0.rds"))
  smr <- my_eval(res_vmrseq$gr_qced)
  # write.csv(smr$true, paste0("code/test/test_vmrseq/true_vmrs_size",r,".csv"))
  # write.csv(smr$detected, paste0("code/test/test_vmrseq/detected_vmrs_size",r,".csv"))
  write.csv(smr$true, paste0("code/test/test_vmrseq/true_vmrs_",N,"cells.csv"))
  write.csv(smr$detected, paste0("code/test/test_vmrseq/detected_vmrs_",N,"cells.csv"))
  # for (penalty in 0:10) {
  #   cat("Penalty =", penalty, "\n")
  #   res_gr <- res_vmrseq$gr_qced
  #   name <- c("vmr_index", "cr_index", "vmr_num_cpg", "optim_pi", "n_iter","loglik_diff")
  #   values(res_gr)[which(res_gr$loglik_diff<penalty), name] <- NA
  #   smr <- my_eval(res_gr)
  #   if (penalty==0) {
  #     write.csv(smr$true, paste0("code/test/test_vmrseq/true_vmrs_",N,"cells.csv"))
  #     write.csv(smr$detected, paste0("code/test/test_vmrseq/detected_vmrs_",N,"cells.csv"))
  #   }
  # }
}
