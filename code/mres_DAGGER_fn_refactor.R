# Functional refactoring of the multi-resolution DAGGER method

library(data.tree)
library(tidyverse)
library(igraph)
#--------------------------
library(GenomicRanges)
library(parallel)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------
# utils Fn.
rej_fn<-function(nodes,lvl,num_rejected,alpha){
  #pick node effective node and leaf number
  ms_d<-ms[nodes]
  ls_d<-ls[nodes]
  p_vals_d<- node_pval[nodes]
  
  
  ### P-value threshold function
  # r is the considered rank
  crit_func <- function(r,alpha){
    alpha * ls_d * (ms_d + r + num_rejected[as.character(lvl-1)] - 1) / l / ms_d
  }
  
  r <- length(p_vals_d)
  # Determine the appropriate threshold value
  while (sum(p_vals_d <= crit_func(r,alpha)) < r){
    r <- r-1 
  }  
  R <- r 
  tmp_rejected<-nodes[which(p_vals_d <= crit_func(R,alpha))]
  
  return(tmp_rejected)
}

detect_inter_cage_cl_fn<-function(feature_coord_tbl,feature_pval_tbl,res_num){
  
  fn_env<-environment()
  
  feature_Grange<-   GRanges(seqnames=feature_coord_tbl$chr,
                          ranges = IRanges(start=feature_coord_tbl$start,
                                           end=feature_coord_tbl$end
                          ))
  
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(GenomicRanges)
    library(dplyr)
    print("node ready")
  })
  clusterExport(cl,c("feature_pval_tbl","feature_Grange","res_num"),envir = fn_env)
  feature_bin_n_l<-parLapply(cl,1:nrow(feature_pval_tbl),function(x){
    cl_Grange<-   GRanges(seqnames=feature_pval_tbl$chr[x],
                          ranges = IRanges(start=as.numeric(feature_pval_tbl$bins[[x]]),
                                           end=as.numeric(feature_pval_tbl$bins[[x]]) + res_num[feature_pval_tbl$res[x]]-1
                          ))
    return(length(unique(queryHits(findOverlaps(cl_Grange,feature_Grange)))))
    
  })
  stopCluster(cl)
  rm(cl)
  feature_pval_tbl<-feature_pval_tbl %>% mutate(feature.bin=unlist(feature_bin_n_l))
  return(feature_pval_tbl)

}

mres_DAGGER_fn<-function(chr_pval_tbl,chr_bpt,alpha_seq){
  
  
}

#--------------------------
feature_coord_file<-"./data/CAGE_tss_coord_HMEC_tbl.Rda"
feature_pval_file<-"./data/pval_tbl/CAGE_tss_HMEC_pval_tbl.Rda"
chromo<-"chr22"

feature_coord_tbl<-get(load(feature_coord_file))
tmp_obj<-names(mget(load(feature_coord_file)))
rm(list=tmp_obj)
rm(tmp_obj)

feature_pval_tbl<-get(load(feature_pval_file))
tmp_obj<-names(mget(load(feature_pval_file)))
rm(list=tmp_obj)
rm(tmp_obj)
detect_inter_cage_cl_fn(feature_coord_tbl %>% filter(!(is.na(start))),feature_pval_tbl,res_num)
