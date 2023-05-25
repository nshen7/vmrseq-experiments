lib_path <- "/home/nshen7/R/rstudio_4_2_0"
if (file.exists(lib_path)) .libPaths(lib_path) else message('Lib path not exist!')
library(tidyverse)
library(data.table)
library(here)
setwd(here())
