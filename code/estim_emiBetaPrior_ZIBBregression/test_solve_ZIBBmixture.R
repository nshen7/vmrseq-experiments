dir <- "data/processed/processed_liu2021/"
file_list <- list.files(dir)
file_name <- file_list[8]
print(file_name)
cells.se <- loadHDF5SummarizedExperiment(dir = paste0(dir, file_name))

sample_n <- 10000; set.seed(2021)
index <- sample(1:nrow(cells.se), sample_n)
sub_dat <- data.frame(
  DataSource = "Liu2021",
  CellClass = info_row$CellClass,
  SubType = info_row$SubType,
  N_cell = ncol(cells.se),
  med_cov = median(values(cells.se)$cell_cov),
  cell_cov = values(cells.se)$cell_cov[index],
  cell_meth = values(cells.se)$cell_meth[index],
  cell_MF = values(cells.se)$cell_MF[index]
) 
sub_dat <- sub_dat %>%
  mutate(Cluster = as.integer(sub_dat$cell_MF >= 0.5)) %>%
  mutate(y_u = as.matrix(data.frame(cell_meth, cell_cov-cell_meth)),
         y_m = as.matrix(data.frame(cell_cov-cell_meth, cell_meth))) 

for (i in 1:10) {
  mod_u <- gamlss(y_u ~ 1, data = sub_dat %>% filter(Cluster==0), family = ZIBB, 
                  trace = F, n.cyc = 100)
  plot(mod_u)
  cat("u", fitted(mod_u, "nu")[1], fitted(mod_u, "mu")[1], unname(fitted(mod_u, "sigma")[1]), "\n")
  mod_m <- gamlss(y_m ~ 1, data = sub_dat %>% filter(Cluster==1), family = ZIBB, 
                  trace = F, n.cyc = 100)
  plot(mod_m)
  cat("m", fitted(mod_m, "nu")[1], fitted(mod_m, "mu")[1], unname(fitted(mod_m, "sigma")[1]), "\n")
  
  c_prior <- unname(table(sub_dat$Cluster) / nrow(sub_dat))
  sub_dat <- sub_dat %>% 
    mutate(p0_u = dZIBB(cell_meth, 
                        nu = fitted(mod_u, "nu")[1], 
                        mu = fitted(mod_u, "mu")[1], 
                        sigma = unname(fitted(mod_u, "sigma")[1]),
                        bd = cell_cov),
           p0_m = dZIBB(cell_cov - cell_meth, 
                        nu = fitted(mod_u, "nu")[1], 
                        mu = fitted(mod_m, "mu")[1], 
                        sigma = unname(fitted(mod_m, "sigma")[1]),
                        bd = cell_cov)) %>%
    mutate(p_u = p0_u*c_prior[1]/(p0_u*c_prior[1] + p0_m*c_prior[2])) %>%
    mutate(Cluster = map_int(p_u, ~ sample(0:1, 1, prob = c(.x, 1-.x))))
  # hist(sub_dat$cell_MF[sub_dat$Cluster==0], xlim = c(0,1), breaks = 100, freq = F)
  # hist(sub_dat$cell_MF[sub_dat$Cluster==1], xlim = c(0,1), breaks = 100, freq = F)
}

## Issue: sites with MF==1 might be assigned to unmethylated cluster.