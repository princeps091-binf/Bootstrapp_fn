library(GenomicRanges)
library(tidyverse)
library(furrr)
library(data.tree)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
get_obj_in_fn<-function(file){
  tmp_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(tmp_tbl)
}
#-----------------------------------------
hub_peak_content_file<-"./data/GM12878_5kb_hub_ancestry_CAGE_peak_content.Rda"
parent_hub_content_file<-"./data/GM12878_5kb_hub_parent_descendance.Rda"
compound_hub_file<-"./data/GM12878_5kb_hub_ancestry.Rda"

hub_peak_content_tbl<-get_obj_in_fn(hub_peak_content_file)
parent_hub_content_tbl<-get_obj_in_fn(parent_hub_content_file)
compound_hub_tbl<-get_obj_in_fn(compound_hub_file)

hub_5kb_cage_content<-hub_peak_content_tbl %>% 
  filter(res=="5kb") %>% 
  dplyr::select(chr,hub,peak.content)

plan(multisession, workers = 3)

parent_hub_content_tbl<-parent_hub_content_tbl %>% 
  group_by(chr,parent.hub) %>% 
  summarise(ch.hub=list(unique(children.hub))) %>% 
  ungroup() %>% 
  mutate(hub.5kb.foot=future_pmap_dbl(list(chr,parent.hub,ch.hub),function(chromo,parent.hub,ch.hub){
    
    parent_peak_content<-hub_peak_content_tbl %>% 
      filter(chr==chromo & hub == parent.hub) %>% 
      dplyr::select(peak.content) %>% unnest(cols=c(peak.content)) %>% distinct %>% unlist
    tmp_ok_hub<-hub_5kb_cage_content %>% filter(chr == chromo) %>% dplyr::select(hub) %>% unlist
    ok_ch<-ch.hub[ch.hub %in% tmp_ok_hub]
    if(length(ok_ch)<1){return(NA)} else{
      child_peak_content<-hub_5kb_cage_content %>% 
        filter(chr==chromo) %>% filter(hub %in% ok_ch) %>% 
        dplyr::select(peak.content) %>% unnest(cols=c(peak.content)) %>% distinct %>% unlist
      return(sum(parent_peak_content %in% child_peak_content)/length(parent_peak_content))
      
    }
    
  })) %>% 
  arrange(desc(hub.5kb.foot))

candidate_hub_tbl<-parent_hub_content_tbl %>% 
  filter(hub.5kb.foot>0.5) %>% 
  mutate(res=map_chr(parent.hub,function(x){
    strsplit(x,split="_")[[1]][1]
  }))


walk(res_num,function(tmp_res){
  message(names(res_num[which(res_num==tmp_res)]))
  tmp_tbl<-candidate_hub_tbl %>% 
    filter(res_num[res] >= tmp_res) 
  
  save(tmp_tbl,file=paste0("./data/candidate_compound_hub/GM12878_",names(res_num[which(res_num==tmp_res)]),"_tss_compound_hub.Rda")) 
  
})

hub_5kb_stub<-compound_hub_tbl %>% 
  filter(is.na(parent.hub)) %>% 
  distinct(chr,hub.5kb) %>% 
  bind_rows(.,compound_hub_tbl %>% 
              filter(!(is.na(parent.hub))) %>% 
              group_by(chr,hub.5kb) %>% 
              summarise(stub=all(grepl("^5kb_",parent.hub))) %>% 
              filter(stub) %>% 
              dplyr::select(-stub) %>% 
              ungroup)
save(hub_5kb_stub,file=paste0("./data/candidate_compound_hub/GM12878_hub_5kb_stub.Rda")) 

