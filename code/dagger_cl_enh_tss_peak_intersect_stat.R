# Examining the coverage of the candidate trx-hubs
library(tidyverse)
library(GenomicRanges)
library(furrr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------
#Utils.Fn,
Get_cl_bin_fn<-function(chr_cl_tbl,res_file,res_num){
    tmp_chromo<- unique(chr_cl_tbl$chr)
    load(paste0(res_file,tmp_chromo,"_spec_res.Rda"))
    
    chr_cl_tbl<-chr_cl_tbl %>% mutate(bin=chr_spec_res$cl_member[node]) 
    

# Output tibble with GRange list-column
  return(chr_cl_tbl)
}

Build_Grange_fn<-function(enh_dagger_tbl,res_file,res_num){
  
  enh_chr_set<-enh_dagger_tbl %>% distinct(chr) %>% unlist
  enh_dagger_cl_bin_tbl<-do.call(bind_rows,lapply(enh_chr_set,function(chromo){
    chr_cl_tbl<-enh_dagger_tbl %>% filter(chr==chromo) %>% select(chr,node)
    Get_cl_bin_fn(chr_cl_tbl,res_file,res_num)
  }))
  
  enh_dagger_tbl<-enh_dagger_tbl %>% left_join(.,enh_dagger_cl_bin_tbl)
  enh_dagger_tbl<-enh_dagger_tbl %>% mutate(end=pmap(list(bin,res),function(bin,res){
    as.numeric(bin) + res_num[res] -1
  }))
  
  plan(multisession, workers = 3)
  enh_dagger_tbl<-enh_dagger_tbl %>% mutate(GRange=future_pmap(list(chr,res,bin,end),function(chr,res,bin,end){
    return(reduce(GRanges(seqnames=chr,
                   ranges = IRanges(start=as.numeric(bin),
                                    end=as.numeric(end)
                   ))))
  }))
  return(enh_dagger_tbl)
}
#--------------------------
# Load the summary tables for candidate trx-hubs
enh_file<-"./data/DAGGER_tbl/HMEC_enh_dagger_tbl.Rda"
tss_file<-"./data/DAGGER_tbl/HMEC_tss_dagger_tbl.Rda"
union_file<-"./data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"

res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"

tss_dagger_tbl<-get(load(tss_file))
tmp_obj<-names(mget(load(tss_file)))
rm(list=tmp_obj)
rm(tmp_obj)
tss_dagger_tbl<-tss_dagger_tbl %>% mutate(ID=paste(chr,node,sep="_"))

enh_dagger_tbl<-get(load(enh_file))
tmp_obj<-names(mget(load(enh_file)))
rm(list=tmp_obj)
rm(tmp_obj)
enh_dagger_tbl<-enh_dagger_tbl %>% mutate(ID=paste(chr,node,sep="_"))

union_dagger_tbl<-get(load(union_file))
tmp_obj<-names(mget(load(union_file)))
rm(list=tmp_obj)
rm(tmp_obj)
union_dagger_tbl<-union_dagger_tbl %>% mutate(ID=paste(chr,node,sep="_"))

enh_dagger_tbl<-enh_dagger_tbl %>% Build_Grange_fn(.,res_file,res_num)

tss_dagger_tbl<-tss_dagger_tbl %>% Build_Grange_fn(.,res_file,res_num)

union_dagger_tbl<-union_dagger_tbl %>% Build_Grange_fn(.,res_file,res_num)

enh_mres_dagger_GRange<-reduce(unlist(GRangesList(enh_dagger_tbl$GRange)))

tss_mres_dagger_GRange<-reduce(unlist(GRangesList(tss_dagger_tbl$GRange)))

union_mres_dagger_GRange<-reduce(unlist(GRangesList(union_dagger_tbl$GRange)))
#--------------------------------------
#Load the corresponding CAGE-tss and CAGE-enh
tss_GRange_file<-"./data/GRanges/CAGE_tss_HMEC_Grange.Rda"
enh_GRange_file<-"./data/GRanges/CAGE_enh_HMEC_Grange.Rda"

tss_GRange<-get(load(tss_GRange_file))
tmp_obj<-names(mget(load(tss_GRange_file)))
rm(list=tmp_obj)
rm(tmp_obj)

enh_GRange<-get(load(enh_GRange_file))
tmp_obj<-names(mget(load(enh_GRange_file)))
rm(list=tmp_obj)
rm(tmp_obj)

length(unique(queryHits(findOverlaps(enh_GRange,enh_mres_dagger_GRange))))/length(enh_GRange)

length(unique(queryHits(findOverlaps(tss_GRange,tss_mres_dagger_GRange))))/length(tss_GRange)

length(unique(queryHits(findOverlaps(enh_GRange,tss_mres_dagger_GRange))))/length(enh_GRange)

length(unique(queryHits(findOverlaps(tss_GRange,enh_mres_dagger_GRange))))/length(tss_GRange)

length(unique(queryHits(findOverlaps(enh_GRange,union_mres_dagger_GRange))))/length(enh_GRange)

length(unique(queryHits(findOverlaps(tss_GRange,union_mres_dagger_GRange))))/length(tss_GRange)

#-------

# Resolution specific examination of the intersection