library(tidyverse)
library(furrr)

#---------------------------------------------------------------------
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
pval_tbl_file<-"./data/pval_tbl/CAGE_union_H1_pval_tbl.Rda"
hub_5kb_ancestry_tbl_file<-"./data/H1_5kb_hub_ancestry.Rda"
cage_peak_Grange_file<-"./data/GRanges/CAGE_union_H1_Grange.Rda"

pval_tbl<-get_obj_in_fn(pval_tbl_file)
compound_hub_tbl<-get_obj_in_fn(hub_5kb_ancestry_tbl_file)
cage_Grange<-get_obj_in_fn(cage_peak_Grange_file)
mcols(cage_Grange)<-tibble(ID=paste("CAGE",1:length(cage_Grange),sep="_"))

# CAGE-peak content for candidate trans-res hubs
hub_peak_content_tbl<- compound_hub_tbl %>%
  distinct(chr,parent.hub) %>% rename(hub=parent.hub) %>% 
  bind_rows(.,  compound_hub_tbl %>% 
              distinct(chr,hub.5kb) %>% rename(hub=hub.5kb)) %>% 
  distinct(chr,hub) %>% 
  filter(!(is.na(hub))) %>% 
  left_join(.,pval_tbl,by=c("hub"="cl","chr"="chr"))
plan(multisession, workers = 3)

hub_peak_content_tbl<-hub_peak_content_tbl %>% 
  mutate(peak.content=future_map(GRange,function(x){
    cage_Grange@elementMetadata$ID[unique(queryHits(findOverlaps(cage_Grange,x)))]
  }))
save(hub_peak_content_tbl,file="./data/H1_5kb_hub_ancestry_CAGE_peak_content.Rda")
