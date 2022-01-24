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

get_data_tbl<-function(file){
  
  data_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(data_tbl)
}

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

tss_dagger_tbl<-get_data_tbl(tss_file)
tss_dagger_tbl<-tss_dagger_tbl %>% mutate(ID=paste(chr,node,sep="_"))

enh_dagger_tbl<-get_data_tbl(enh_file)
enh_dagger_tbl<-enh_dagger_tbl %>% mutate(ID=paste(chr,node,sep="_"))

union_dagger_tbl<-get_data_tbl(union_file)
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

tss_GRange<-get_data_tbl(tss_GRange_file)

enh_GRange<-get_data_tbl(enh_GRange_file)


length(unique(queryHits(findOverlaps(enh_GRange,enh_mres_dagger_GRange))))/length(enh_GRange)

length(unique(queryHits(findOverlaps(tss_GRange,tss_mres_dagger_GRange))))/length(tss_GRange)

length(unique(queryHits(findOverlaps(enh_GRange,tss_mres_dagger_GRange))))/length(enh_GRange)

length(unique(queryHits(findOverlaps(tss_GRange,enh_mres_dagger_GRange))))/length(tss_GRange)

length(unique(queryHits(findOverlaps(enh_GRange,union_mres_dagger_GRange))))/length(enh_GRange)

length(unique(queryHits(findOverlaps(tss_GRange,union_mres_dagger_GRange))))/length(tss_GRange)

#-------

# Resolution specific examination of the intersection

dagger_res_set<-unique(c(enh_dagger_tbl$res,tss_dagger_tbl$res,union_dagger_tbl$res))
#enhancer
enh_over_l<-lapply(dagger_res_set,function(tmp_res){
  message(tmp_res)
  tmp_enh_dagger_tbl<-enh_dagger_tbl %>% filter(res==tmp_res)
  tmp_tss_dagger_tbl<-tss_dagger_tbl %>% filter(res==tmp_res)
  tmp_union_dagger_tbl<-union_dagger_tbl %>% filter(res==tmp_res)
  
  
  tmp_enh_dagger_GRange<-reduce(unlist(GRangesList(tmp_enh_dagger_tbl$GRange)))
  
  tmp_tss_dagger_GRange<-reduce(unlist(GRangesList(tmp_tss_dagger_tbl$GRange)))
  
  tmp_union_dagger_GRange<-reduce(unlist(GRangesList(tmp_union_dagger_tbl$GRange)))
  
  enh_over<-length(unique(queryHits(findOverlaps(enh_GRange,tmp_enh_dagger_GRange))))
  tss_over<-length(unique(queryHits(findOverlaps(enh_GRange,tmp_tss_dagger_GRange))))
  uni_over<-length(unique(queryHits(findOverlaps(enh_GRange,tmp_union_dagger_GRange))))
  n_vec<-c(enh_over,length(enh_GRange)-enh_over,tss_over,length(enh_GRange)-tss_over,uni_over,length(enh_GRange)-uni_over)
  return(tibble(res=tmp_res,hub=rep(c("enh","tss","union"),each=2),io=rep(c("in","out"),3),n=n_vec))
  
})

enh_over_tbl<-do.call(bind_rows,enh_over_l)

mres_enh_over<-length(unique(queryHits(findOverlaps(enh_GRange,enh_mres_dagger_GRange))))
mres_tss_over<-length(unique(queryHits(findOverlaps(enh_GRange,tss_mres_dagger_GRange))))
mres_uni_over<-length(unique(queryHits(findOverlaps(enh_GRange,union_mres_dagger_GRange))))
n_vec<-c(mres_enh_over,length(enh_GRange)-mres_enh_over,mres_tss_over,length(enh_GRange)-mres_tss_over,mres_uni_over,length(enh_GRange)-mres_uni_over)
enh_over_tbl<-enh_over_tbl %>% bind_rows(.,tibble(tibble(res="mres",hub=rep(c("enh","tss","union"),each=2),io=rep(c("in","out"),3),n=n_vec)))
enh_over_tbl %>%
  mutate(res=fct_relevel(res,c("mres",rev(names(sort(res_num)))))) %>% 
  ggplot(.,aes(hub,n,fill=io))+geom_bar(stat="identity")+
  facet_wrap(res~.)+
  scale_fill_brewer(palette="Set1")

ggsave("~/Documents/multires_bhicect/weeklies/weekly48/img/enh_dagger_coverage.svg")
#tss

tss_over_l<-lapply(dagger_res_set,function(tmp_res){
  message(tmp_res)
  tmp_enh_dagger_tbl<-enh_dagger_tbl %>% filter(res==tmp_res)
  tmp_tss_dagger_tbl<-tss_dagger_tbl %>% filter(res==tmp_res)
  tmp_union_dagger_tbl<-union_dagger_tbl %>% filter(res==tmp_res)
  
  
  tmp_enh_dagger_GRange<-reduce(unlist(GRangesList(tmp_enh_dagger_tbl$GRange)))
  
  tmp_tss_dagger_GRange<-reduce(unlist(GRangesList(tmp_tss_dagger_tbl$GRange)))
  
  tmp_union_dagger_GRange<-reduce(unlist(GRangesList(tmp_union_dagger_tbl$GRange)))
  
  enh_over<-length(unique(queryHits(findOverlaps(tss_GRange,tmp_enh_dagger_GRange))))
  tss_over<-length(unique(queryHits(findOverlaps(tss_GRange,tmp_tss_dagger_GRange))))
  uni_over<-length(unique(queryHits(findOverlaps(tss_GRange,tmp_union_dagger_GRange))))
  n_vec<-c(enh_over,length(tss_GRange)-enh_over,tss_over,length(tss_GRange)-tss_over,uni_over,length(tss_GRange)-uni_over)
  return(tibble(res=tmp_res,hub=rep(c("enh","tss","union"),each=2),io=rep(c("in","out"),3),n=n_vec))
  
})

tss_over_tbl<-do.call(bind_rows,tss_over_l)

mres_enh_over<-length(unique(queryHits(findOverlaps(tss_GRange,enh_mres_dagger_GRange))))
mres_tss_over<-length(unique(queryHits(findOverlaps(tss_GRange,tss_mres_dagger_GRange))))
mres_uni_over<-length(unique(queryHits(findOverlaps(tss_GRange,union_mres_dagger_GRange))))
n_vec<-c(mres_enh_over,length(tss_GRange)-mres_enh_over,mres_tss_over,length(tss_GRange)-mres_tss_over,mres_uni_over,length(tss_GRange)-mres_uni_over)
tss_over_tbl<-tss_over_tbl %>% bind_rows(.,tibble(tibble(res="mres",hub=rep(c("enh","tss","union"),each=2),io=rep(c("in","out"),3),n=n_vec)))
tss_over_tbl %>%
  mutate(res=fct_relevel(res,c("mres",rev(names(sort(res_num)))))) %>% 
  ggplot(.,aes(hub,n,fill=io))+geom_bar(stat="identity")+
  facet_wrap(res~.)+
  scale_fill_brewer(palette="Set1")
ggsave("~/Documents/multires_bhicect/weeklies/weekly48/img/tss_dagger_coverage.svg")

