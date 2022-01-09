# Produce Coordinate tables for TADs
renv::install("tidyverse")
renv::install("bioc::GenomicRanges")
renv::install("furrr")

library(tidyverse)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------------------------------------------------
## Import table with Feature coordinate
feature_file<-"~/Documents/multires_bhicect/data/HMEC/GSE63525_HMEC_Arrowhead_domainlist.txt.gz"
out_folder<-"./data/TAD_GRange/"
feature_tbl<-read_delim(feature_file,
                                   delim = "\t", escape_double = FALSE,
                                   col_names = T, trim_ws = TRUE)
feature_tbl<-feature_tbl%>%mutate(chr1=paste0("chr",chr1))
plan(multisession, workers = 3)

feature_tbl<-feature_tbl %>% 
    mutate(GRange=future_pmap(list(chr1,x1,x2),function(chr1,x1,x2){
      return(GRanges(seqnames=chr1,
                              ranges = IRanges(start=x1,
                                               end=x2
                              )))
  
})) %>% 
  dplyr::select(chr1,x1,x2,GRange)%>% 
  tidyr::unite(ID,chr1,x1,x2,sep="_",remove=F) %>% 
  dplyr::select(ID,chr1,GRange) %>% dplyr::rename(chr=chr1)
lapply(unique(feature_tbl$chr),function(chromo){
  chr_tbl<-feature_tbl %>% filter(chr==chromo)
  save(chr_tbl,file=paste0(out_folder,chromo,"_TAD.Rda"))
  
  
})
