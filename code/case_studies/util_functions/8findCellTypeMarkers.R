source("code/SETPATHS.R")

findCellTypeMarkerIndices <- function(target_type, target_type_n_cell, ct.se, min_diff_hyper, min_diff_hypo, min_avail_cell_prop = 0.2) {
  
  b_tg <- assays(ct.se)$reginal_methyl[,  which(colnames(ct.se) == target_type)]
  B_bg <- assays(ct.se)$reginal_methyl[, -which(colnames(ct.se) == target_type)]
  
  diff_hyper <- b_tg - rowMaxs(as.matrix(B_bg), na.rm = T)
  diff_hypo  <- b_tg - rowMins(as.matrix(B_bg), na.rm = T)
  
  cell_avail <- assays(ct.se)$n_cell[, which(colnames(ct.se) == target_type)]
  incl_min_avail_cell_prop <- cell_avail > min_avail_cell_prop * target_type_n_cell
  idx_hyper <- which(diff_hyper >  abs(min_diff_hyper))
  idx_hypo  <- which(diff_hypo  < -abs(min_diff_hypo))
  
  return(list(idx_hyper = idx_hyper, idx_hypo = idx_hypo))
}