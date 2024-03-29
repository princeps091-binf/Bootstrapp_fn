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
pval_tbl_file<-"./data/pval_tbl/CAGE_union_HMEC_pval_tbl.Rda"
cage_peak_Grange_file<-"./data/GRanges/CAGE_union_HMEC_Grange.Rda"
hub_5kb_ancestry_file<-"./data/HMEC_5kb_hub_ancestry.Rda"
hub_ancestry_peak_content_file<-"./data/HMEC_5kb_hub_ancestry_CAGE_peak_content.Rda"
cage_Grange<-get_obj_in_fn(cage_peak_Grange_file)
mcols(cage_Grange)<-tibble(ID=paste("CAGE",1:length(cage_Grange),sep="_"))

pval_tbl<-get_obj_in_fn(pval_tbl_file)
dagger_mres_hub_tbl<-get_obj_in_fn(dagger_union_file)
compound_hub_tbl<-get_obj_in_fn(hub_5kb_ancestry_file)
hub_peak_content_tbl<-get_obj_in_fn(hub_ancestry_peak_content_file)
#---------------------------------------------------------------------

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

walk(res_num,function(tmp_res){
  
  hub_5kb_gene<-candidate_hub_tbl %>% 
    filter(res_num[res] >= tmp_res) %>% 
    unnest(cols=c(ch.hub)) %>% 
    mutate(ch.res=map_chr(ch.hub,function(x){
      strsplit(x,split="_")[[1]][1]
    })) %>% 
    filter(ch.res=="5kb") %>% 
    distinct(chr,ch.hub)  %>% 
    inner_join(.,hub_5kb_cage_content,by=c("chr"="chr","ch.hub"="hub")) %>% 
    unnest(cols=c(peak.content)) %>% distinct(peak.content) %>% unlist
  
  
  full_gene<-candidate_hub_tbl %>% 
    filter(res_num[res] >= tmp_res) %>% 
    unnest(cols=c(ch.hub)) %>% 
    inner_join(.,hub_peak_content_tbl,by=c("chr"="chr","parent.hub"="hub")) %>% 
    unnest(cols=c(peak.content)) %>% distinct(peak.content) %>% unlist
  message(sum(full_gene %in% hub_5kb_gene)/length(unique(c(full_gene))))
})

ori_5kb_hub_content<-hub_5kb_cage_content %>% 
  unnest(cols=c(peak.content)) %>% distinct(peak.content) %>% unlist


walk(res_num,function(tmp_res){
  
  hub_5kb_gene<-candidate_hub_tbl %>% 
    filter(res_num[res] >= tmp_res) %>% 
    unnest(cols=c(ch.hub)) %>% 
    mutate(ch.res=map_chr(ch.hub,function(x){
      strsplit(x,split="_")[[1]][1]
    })) %>% 
    filter(ch.res=="5kb") %>% 
    distinct(chr,ch.hub)  %>% 
    inner_join(.,hub_5kb_cage_content,by=c("chr"="chr","ch.hub"="hub")) %>% 
    unnest(cols=c(peak.content)) %>% distinct(peak.content) %>% unlist
  
  message(sum(hub_5kb_gene %in% ori_5kb_hub_content)/length(unique(c(ori_5kb_hub_content))))
})

walk(res_num,function(tmp_res){
  
tmp_peak<-candidate_hub_tbl %>% 
  filter(res_num[res] >= tmp_res) %>% 
  left_join(.,hub_peak_content_tbl %>% dplyr::select(chr,hub,peak.content),by=c("chr"="chr","parent.hub"="hub")) %>% 
  unnest(cols=c(peak.content)) %>% 
  distinct(peak.content) %>% unlist
message(sum(tmp_peak %in% cage_Grange@elementMetadata$ID)/length(cage_Grange@elementMetadata$ID))
})


hub_5kb_cage<-dagger_mres_hub_tbl%>% 
  filter(res=="5kb") %>% 
  inner_join(.,hub_peak_content_tbl %>% dplyr::select(chr,hub,peak.content),by=c("chr"="chr","node"="hub")) %>% 
  unnest(cols=c(peak.content)) %>% 
  distinct(peak.content) %>% unlist

sum(hub_5kb_cage %in% cage_Grange@elementMetadata$ID)/length(cage_Grange@elementMetadata$ID)

walk(res_num,function(tmp_res){
  message(names(res_num[which(res_num==tmp_res)]))
  tmp_tbl<-candidate_hub_tbl %>% 
    filter(res_num[res] >= tmp_res) %>% 
    left_join(.,hub_peak_content_tbl %>% dplyr::select(chr,hub,peak.content),by=c("chr"="chr","parent.hub"="hub"))
  
  save(tmp_tbl,file=paste0("./data/candidate_compound_hub/HMEC_",names(res_num[which(res_num==tmp_res)]),"_tss_compound_hub.Rda")) 
    
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
save(hub_5kb_stub,file=paste0("./data/candidate_compound_hub/HMEC_hub_5kb_stub.Rda")) 

hub_5kb_stub_cage<-hub_5kb_stub %>% 
  inner_join(.,hub_peak_content_tbl %>% dplyr::select(chr,hub,peak.content),by=c("chr"="chr","hub.5kb"="hub")) %>% 
  unnest(cols=c(peak.content)) %>% 
  distinct(peak.content) %>% unlist

sum(hub_5kb_stub_cage %in% cage_Grange@elementMetadata$ID)/length(cage_Grange@elementMetadata$ID)
