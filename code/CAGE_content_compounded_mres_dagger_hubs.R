# Examine the gene content of mres-compound DAGGER hubs compared to single res DAGGER hubs
library(tidyverse)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-------------------------------
get_obj_in_fn<-function(file){
  tmp_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(tmp_tbl)
  }
##input hic data of choice

hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat <- read_delim(paste0(dat_file,cl_res,"/",chromo,".txt"), 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.nan(X3)))%>%filter(X1!=X2))
}

Build_GRange_fn<-function(bin_set,res,chr,res_num){
  
  return(GRanges(seqnames=chr,
                 ranges = IRanges(start=as.numeric(bin_set),
                                  end=as.numeric(bin_set)+res_num[res]-1
                 )))
  
}
#-------------------------------
tss_Grange_file<-"./data/GRanges/CAGE_tss_HMEC_Grange.Rda"
compound_mres_hub_file<-"./data/DAGGER_tbl/HMEC_mres_coumpund_hub.Rda"
pval_tbl_file<-"./data/pval_tbl/CAGE_union_HMEC_pval_tbl.Rda"
dagger_union_file<-"./data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"

tss_Grange<-get_obj_in_fn(tss_Grange_file)
compound_mres_hub_tbl<-get_obj_in_fn(compound_mres_hub_file)
pval_tbl<-get_obj_in_fn(pval_tbl_file)
dagger_mres_hub_tbl<-get_obj_in_fn(dagger_union_file)

top_hub_tbl<-compound_mres_hub_tbl %>% 
  distinct(chr,top_hub,top_res) %>% 
  filter(!(is.na(top_hub)))

top_hub_tbl<-top_hub_tbl %>% 
  left_join(.,pval_tbl %>% 
              dplyr::select(chr,cl,GRange,emp.pval),by=c("chr"='chr','top_hub'="cl"))

top_hub_GRange<-GenomicRanges::reduce(top_hub_tbl %>% 
  filter(top_res %in% c("1Mb","500kb")) %>% dplyr::select(GRange) %>% unlist %>% GRangesList %>% unlist)

lores_hub_GRange<-GenomicRanges::reduce(dagger_mres_hub_tbl %>% 
                                        filter(res %in% c("1Mb","500kb")) %>% 
                                        left_join(.,pval_tbl %>% 
                                                    dplyr::select(chr,cl,GRange,emp.pval),by=c("chr"='chr','node'="cl")) %>% 
                                        dplyr::select(GRange) %>% unlist %>% GRangesList %>% unlist)

hires_hub_GRange<-GenomicRanges::reduce(dagger_mres_hub_tbl %>% 
                                          filter(res %in% c("5kb")) %>% 
                                          left_join(.,pval_tbl %>% 
                                                      dplyr::select(chr,cl,GRange,emp.pval),by=c("chr"='chr','node'="cl")) %>% 
                                          dplyr::select(GRange) %>% unlist %>% GRangesList %>% unlist)


tss_Grange[unique(subjectHits(findOverlaps(top_hub_GRange,tss_Grange)))]

tss_Grange[unique(subjectHits(findOverlaps(lores_hub_GRange,tss_Grange)))]

tss_Grange[unique(subjectHits(findOverlaps(hires_hub_GRange,tss_Grange)))]
