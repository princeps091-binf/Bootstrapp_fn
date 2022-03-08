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
  out_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#-----------------------------------------
compound_hub_10kb_file<-"./data/candidate_compound_hub/GM12878_10kb_tss_compound_hub.Rda"
compound_hub_5kb_file<-"./data/candidate_compound_hub/GM12878_5kb_tss_compound_hub.Rda"
compound_hub_file<-"./data/GM12878_5kb_hub_ancestry.Rda"

compound_hub_tbl<-get_obj_in_fn(compound_hub_file)
compound_hub_5kb_tbl<-get_obj_in_fn(compound_hub_5kb_file)
compound_hub_10kb_tbl<-get_obj_in_fn(compound_hub_10kb_file)


hub_5kb_stub<-compound_hub_tbl %>% 
  group_by(chr,hub.5kb) %>% 
  summarise(is.5kb.stub=all(grepl("^5kb_",parent.hub) | is.na(parent.hub))) %>% 
  filter(is.5kb.stub)

non_5kb_stub<-compound_hub_tbl %>% 
  group_by(chr,hub.5kb) %>% 
  summarise(is.5kb.stub=all(grepl("^5kb_",parent.hub) | is.na(parent.hub))) %>% 
  filter(!(is.5kb.stub)) %>% 
  left_join(.,compound_hub_tbl) 

non_5kb_stub_tbl<-non_5kb_stub %>% 
  distinct(chr,parent.hub) %>% 
  dplyr::rename(hub=parent.hub) %>% 
  bind_rows(.,non_5kb_stub %>% 
              distinct(chr,hub.5kb) %>% 
              dplyr::rename(hub=hub.5kb)) %>% 
  distinct(chr,hub)

rm(non_5kb_stub)
compound_hub_5kb_tbl %>% 
  inner_join(.,hub_5kb_stub,by=c("chr"="chr","parent.hub"="hub.5kb")) %>% 
  arrange(hub.5kb.foot)
  
