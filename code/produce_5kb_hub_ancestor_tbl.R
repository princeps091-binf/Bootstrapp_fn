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

#-------------------------------
pval_tbl_file<-"./data/pval_tbl/CAGE_union_H1_pval_tbl.Rda"
dagger_union_file<-"./data/DAGGER_tbl/H1_union_dagger_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/H1/Dekker/spec_res/"

dagger_mres_hub_tbl<-get_obj_in_fn(dagger_union_file)

tmp_res<-"5kb"
chr_set<-unique(dagger_mres_hub_tbl$chr)
compound_hub_l<-vector("list",length(chr_set))
names(compound_hub_l)<-chr_set


for(chromo in unique(dagger_mres_hub_tbl$chr)){
  
  compound_hub_l[[chromo]]<-get_tophub_fn(dagger_mres_hub_tbl,spec_res_file,chromo,tmp_res)
  
}
compound_hub_tbl<-do.call(bind_rows,compound_hub_l)

save(compound_hub_tbl,file="./data/H1_5kb_hub_ancestry.Rda")
compound_hub_tbl<-get_obj_in_fn("./data/HMEC_5kb_hub_ancestry.Rda")
