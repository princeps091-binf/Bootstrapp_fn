# Comparing intersection between different forms of enrichment
library(tidyverse)
library(UpSetR)
library(svglite)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------

enh_file<-"./data/DAGGER_tbl/HMEC_enh_dagger_tbl.Rda"
tss_file<-"./data/DAGGER_tbl/HMEC_tss_dagger_tbl.Rda"
union_file<-"./data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"


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

listInput<-list(enh=enh_dagger_tbl$ID,tss=tss_dagger_tbl$ID,union=union_dagger_tbl$ID)
  
svglite(filename = paste0("./img/dagger_cl_upset_mres.svg"))
upset(fromList(listInput), order.by = "freq")
dev.off()

tmp_res_set<-unique(union_dagger_tbl$res)  
walk(tmp_res_set,function(tmp_res){
  tmp_enh_tbl<-enh_dagger_tbl %>% filter(res==tmp_res)
  tmp_tss_tbl<-tss_dagger_tbl %>% filter(res==tmp_res)
  tmp_uni_tbl<-union_dagger_tbl %>% filter(res==tmp_res)
  listInput<-list(enh=tmp_enh_tbl$ID,tss=tmp_tss_tbl$ID,union=tmp_uni_tbl$ID)
  
  svglite(filename = paste0("./img/dagger_cl_upset_",tmp_res,".svg"))
  print({
  upset(fromList(listInput), order.by = "freq")
  })
  dev.off()
  
})

listInput<-list(enh=enh_dagger_tbl$ID,tss=tss_dagger_tbl$ID)

svglite(filename = paste0("./img/dagger_cl_upset_enh_tss_mres.svg"))
upset(fromList(listInput), order.by = "freq")
dev.off()

tmp_res_set<-unique(union_dagger_tbl$res)  
walk(tmp_res_set,function(tmp_res){
  tmp_enh_tbl<-enh_dagger_tbl %>% filter(res==tmp_res)
  tmp_tss_tbl<-tss_dagger_tbl %>% filter(res==tmp_res)
  listInput<-list(enh=tmp_enh_tbl$ID,tss=tmp_tss_tbl$ID)
  
  svglite(filename = paste0("./img/dagger_cl_upset_enh_tss_",tmp_res,".svg"))
  print({
    upset(fromList(listInput), order.by = "freq")
  })
  dev.off()
  
})
