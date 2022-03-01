# Trans-resolution trx-enrichment hubs using CAGE-peak content
library(GenomicRanges)
library(tidyverse)
library(furrr)
library(data.tree)
library(igraph)
library(furrr)
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

Build_GRange_fn<-function(bin_set,res,chr,res_num){
  
  return(GRanges(seqnames=chr,
                 ranges = IRanges(start=as.numeric(bin_set),
                                  end=as.numeric(bin_set)+res_num[res]-1
                 )))
  
}

get_tophub_fn<-function(dagger_hub_tbl,spec_res_file,chromo,tmp_res){
  message(chromo)
  # Get highest resolution seed hubs
  tmp_hub_set<-dagger_hub_tbl %>% filter(chr==chromo & res==tmp_res)
  chr_hubs<-dagger_hub_tbl %>% filter(chr==chromo ) %>% distinct(node) %>% unlist
  # Load the corresponding BPT
  base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  # Collect node ancestors
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  # Collect node levels
  node_lvl<-chr_bpt$Get("level")
  # for each high-resolution hubs, collect the parent hubs
  compound_hubs<-do.call(bind_rows,lapply(tmp_hub_set$node,function(x){
    if(sum(node_ancestor[[x]] %in% chr_hubs)>0){
      tmp_vec<-node_ancestor[[x]]
      return(tibble(chr=chromo,hub.5kb=x,parent.hub=tmp_vec[which(tmp_vec %in% chr_hubs)],parent.hub.lvl=node_lvl[tmp_vec[which(tmp_vec %in% chr_hubs)]]))
    } else{
      return(tibble(chr=chromo,hub.5kb=x,parent.hub=NA,parent.hub.lvl=NA))
    }
    
  }))
  
  return(compound_hubs)
  
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
pval_tbl_file<-"./data/pval_tbl/CAGE_union_HMEC_pval_tbl.Rda"
dagger_union_file<-"./data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
cage_peak_Grange_file<-"./data/GRanges/CAGE_union_HMEC_Grange.Rda"

pval_tbl<-get_obj_in_fn(pval_tbl_file)
dagger_mres_hub_tbl<-get_obj_in_fn(dagger_union_file)
cage_Grange<-get_obj_in_fn(cage_peak_Grange_file)
mcols(cage_Grange)<-tibble(ID=paste("CAGE",1:length(cage_Grange),sep="_"))

dagger_mres_hub_tbl<-dagger_mres_hub_tbl %>% 
  left_join(.,pval_tbl,by=c("node"="cl","chr"="chr","res"="res","emp.pval"="emp.pval"))

tmp_res<-"5kb"
chr_set<-unique(dagger_mres_hub_tbl$chr)
compound_hub_l<-vector("list",length(chr_set))
names(compound_hub_l)<-chr_set


for(chromo in unique(dagger_mres_hub_tbl$chr)){
  
  compound_hub_l[[chromo]]<-get_tophub_fn(dagger_mres_hub_tbl,spec_res_file,chromo,tmp_res)
  
}
compound_hub_tbl<-do.call(bind_rows,compound_hub_l)

parent_hub_tbl<-compound_hub_tbl %>%
  filter(!(is.na(parent.hub))) %>% 
  distinct(chr,parent.hub,parent.hub.lvl)

parent_hub_tbl<-parent_hub_tbl %>% 
  left_join(.,pval_tbl,by=c("parent.hub"="cl","chr"="chr"))
#---------------------------------------------------------------------
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
#---------------------------------------------------------------------
# Compute the CAGE-peak redundancy of trans-res hubs

## Collect for each parent-hub the nested hub content
chr_set<-unique(parent_hub_tbl$chr)
parent_hub_content_l<-vector("list",length(chr_set))
names(parent_hub_content_l)<-chr_set

for(chromo in chr_set){
  
  parent_hub_content_l[[chromo]]<-get_children_hub(parent_hub_tbl,dagger_mres_hub_tbl,chromo,spec_res_file)
  
}

parent_hub_content_tbl<-do.call(bind_rows,parent_hub_content_l)


tmp_tbl<-parent_hub_content_tbl  %>%
  group_by(chr,parent.hub) %>% 
  summarise(ch.hub=list(unique(children.hub))) %>% 
  ungroup() 
tmp_tbl<-tmp_tbl %>% 
  mutate(parent.res=map_chr(parent.hub,function(x){
    return(strsplit(x,split="_")[[1]][1])
  }))
min_res<-"5kb"
min_res_tbl<-tmp_tbl %>% 
  filter(parent.res==min_res)
lower_res_tbl<-tmp_tbl %>% 
  filter(parent.res!=min_res)
## Loop lower-resolution in stepwise manner to aggregate hubs composed mostly by higher-res hubs
trans_res_hub_l<-vector('list',length(res_num))
names(trans_res_hub_l)<-names(sort(res_num))
trans_res_hub_l[[min_res]]<-min_res_tbl

tmp_ok_hub_tbl<-do.call(bind_rows,lapply(trans_res_hub_l,function(x){
  if(!(is.null(x))){
    return(x %>% 
             unnest(ch.hub) %>% 
             distinct(chr,parent.hub) %>% rename(hub=parent.hub) %>% 
             bind_rows(.,  x %>% 
                         unnest(ch.hub) %>% 
                         distinct(chr,ch.hub) %>% rename(hub=ch.hub)) %>% 
             distinct(chr,hub))    
  } else{
    return(tibble(chr=NA,hub=NA))
  }
  
  
})) %>% filter(!(is.na(hub)))


for (tmp_res in names(sort(res_num[which(res_num > res_num[min_res])]))){
  message(tmp_res)
  tmp_res_tbl<-lower_res_tbl %>%
    filter(parent.res == tmp_res) %>% 
    mutate(hub.foot=pmap_dbl(list(chr,parent.hub,ch.hub),function(chromo,parent.hub,ch.hub){
      parent_peak_content<-hub_peak_content_tbl %>% 
        filter(chr==chromo & hub == parent.hub) %>% 
        dplyr::select(peak.content) %>% unnest(cols=c(peak.content)) %>% distinct %>% unlist
      tmp_ok_hub<-tmp_ok_hub_tbl %>% filter(chr == chromo) %>% dplyr::select(hub) %>% unlist
      ok_ch<-ch.hub[ch.hub %in% tmp_ok_hub]
      if(length(ok_ch)<1){return(NA)} else{
        child_peak_content<-hub_peak_content_tbl %>% 
          filter(chr==chromo) %>% filter(hub %in% ok_ch) %>% 
          dplyr::select(peak.content) %>% unnest(cols=c(peak.content)) %>% distinct %>% unlist
        return(sum(parent_peak_content %in% child_peak_content)/length(parent_peak_content))
        
      }
    }))
  trans_res_hub_l[[tmp_res]]<-tmp_res_tbl %>% filter(hub.foot > 0.5)
  # Collect ALL higher-resolution valid hubs
  tmp_ok_hub_tbl<-do.call(bind_rows,lapply(trans_res_hub_l,function(x){
    if(!(is.null(x))){
      return(x %>% 
               unnest(ch.hub) %>% 
               distinct(chr,parent.hub) %>% rename(hub=parent.hub) %>% 
               bind_rows(.,  x %>% 
                           unnest(ch.hub) %>% 
                           distinct(chr,ch.hub) %>% rename(hub=ch.hub)) %>% 
               distinct(chr,hub))    
    } else{
      return(tibble(chr=NA,hub=NA))
    }
    
    
  })) %>% filter(!(is.na(hub)))
  
}
