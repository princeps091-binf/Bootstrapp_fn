library(GenomicRanges)
library(tidyverse)
library(furrr)
library(data.tree)
library(igraph)
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
get_children_hub<-function(parent_hub_tbl,dagger_mres_hub_tbl,chromo,spec_res_file){
  
  message(chromo)
  # Collect parent hub of interest
  tmp_hub_set<-parent_hub_tbl %>% filter(chr==chromo)
  chr_hubs<-dagger_mres_hub_tbl %>% filter(chr==chromo ) %>% distinct(node) %>% unlist
  # Load the corresponding BPT
  base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  # Collect all the children hubs for the considered parent hubs
  compound_hubs<-do.call(bind_rows,lapply(tmp_hub_set$parent.hub,function(x){
    tmp_ch<-names(which(unlist(lapply(node_ancestor[chr_hubs],function(y) x %in% y))))
    return(tibble(chr=chromo,parent.hub=x,children.hub=tmp_ch))
  }))
  
  return(compound_hubs)
  
}
#-------------------------------
dagger_union_file<-"./data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"

dagger_mres_hub_tbl<-get_obj_in_fn(dagger_union_file)
compound_hub_tbl<-get_obj_in_fn("./data/HMEC_5kb_hub_ancestry.Rda")
hub_peak_content_tbl<-get_obj_in_fn("./data/HMEC_5kb_hub_ancestry_CAGE_peak_content.Rda")
#---------------------------------------------------------------------
compound_hub_tbl %>% 
  filter(!(is.na(parent.hub))) %>% 
  group_by(chr,hub.5kb) %>% 
  summarise(lower.res=length(base::grep("5kb",parent.hub,invert=T))>0) %>% 
  filter(!(lower.res))

parent_hub_tbl<-compound_hub_tbl %>%
  filter(!(is.na(parent.hub))) %>% 
  distinct(chr,parent.hub,parent.hub.lvl)

parent_hub_tbl<-parent_hub_tbl %>% 
  left_join(.,pval_tbl,by=c("parent.hub"="cl","chr"="chr"))

# Compute the CAGE-peak redundancy of trans-res hubs

## Collect for each parent-hub the nested hub content
chr_set<-unique(parent_hub_tbl$chr)
parent_hub_content_l<-vector("list",length(chr_set))
names(parent_hub_content_l)<-chr_set

for(chromo in chr_set){
  
  parent_hub_content_l[[chromo]]<-get_children_hub(parent_hub_tbl,dagger_mres_hub_tbl,chromo,spec_res_file)
  
}

parent_hub_content_tbl<-do.call(bind_rows,parent_hub_content_l)

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

candidate_hub_tbl %>% 
  ggplot(.,aes(hub.5kb.foot))+
  geom_density()+
  facet_wrap(res~.,scales="free")
lapply(res_num,function(tmp_res){
  
  hub_5kb_gene<-candidate_hub_tbl %>% 
    filter(res_num[res] >= res_num["10kb"]) %>% 
    unnest(cols=c(ch.hub)) %>% 
    mutate(ch.res=map_chr(ch.hub,function(x){
      strsplit(x,split="_")[[1]][1]
    })) %>% 
    filter(ch.res=="5kb") %>% 
    distinct(chr,ch.hub)  %>% 
    inner_join(.,hub_5kb_cage_content,by=c("chr"="chr","ch.hub"="hub")) %>% 
    unnest(cols=c(peak.content)) %>% distinct(peak.content) %>% unlist
  
  
  full_gene<-candidate_hub_tbl %>% 
    filter(res_num[res] >= res_num["10kb"]) %>% 
    unnest(cols=c(ch.hub)) %>% 
    inner_join(.,hub_peak_content_tbl,by=c("chr"="chr","parent.hub"="hub")) %>% 
    unnest(cols=c(peak.content)) %>% distinct(peak.content) %>% unlist
  sum(full_gene %in% hub_5kb_gene)/length(full_gene)
})
candidate_hub_tbl %>% 
  filter(res_num[res] >= res_num["50kb"]) %>% 
  unnest(cols=c(ch.hub)) %>% 
  mutate(ch.res=map_chr(ch.hub,function(x){
    strsplit(x,split="_")[[1]][1]
  })) %>% 
  filter(ch.res=="5kb") %>% 
  distinct(chr,ch.hub)  %>% 
  inner_join(.,hub_5kb_cage_content,by=c("chr"="chr","ch.hub"="hub")) %>% 
  unnest(cols=c(peak.content)) %>% distinct(peak.content)


candidate_hub_tbl %>% 
  filter(res_num[res] >= res_num["50kb"]) %>% 
  unnest(cols=c(ch.hub)) %>% 
  inner_join(.,hub_peak_content_tbl,by=c("chr"="chr","parent.hub"="hub")) %>% 
  unnest(cols=c(peak.content)) %>% distinct(peak.content)

hub_5kb_cage_content %>% 
  unnest(cols=c(peak.content)) %>% distinct(peak.content)
