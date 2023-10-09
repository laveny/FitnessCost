rm(list = ls())
library(parallel)
setwd("./split.pacbio.bar/")

files <- list.files()

out_path = "./output/"

in_path = "./split.pacbio.bar/"

log_path = "./log/"

mclapply(files,function(x){
  
  system(paste('nohup time python create_all_sim_bar_illu.py ',
               x,' ',in_path,' ',out_path,' > ',log_path,x,'.log 2>&1 &',sep=''))
  
},mc.cores = 50)
