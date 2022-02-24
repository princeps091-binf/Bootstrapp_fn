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
input_data_fn<-function(tmp_file){
  tmp_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(tmp_tbl) 
}

build_Grange<-function(tmp_tbl){
  
  cage_chr_Grange<-GRanges(seqnames=tmp_tbl$chr,
                           ranges = IRanges(start=as.numeric(tmp_tbl$start),
                                            end=as.numeric(tmp_tbl$end)
                           ))
  return(cage_chr_Grange)
  
}

get_tophub_fn<-function(dagger_hub_tbl,spec_res_file,chromo,tmp_res){
  message(chromo)
  tmp_hub_set<-dagger_hub_tbl %>% filter(chr==chromo & res==tmp_res)
  chr_hubs<-dagger_hub_tbl %>% filter(chr==chromo ) %>% distinct(node) %>% unlist
  
  base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  compound_hubs<-unlist(lapply(node_ancestor[tmp_hub_set$node],function(x){
    if(sum(x %in% chr_hubs)>0){
      return(x[max(which(x %in% chr_hubs))])
    } else{
      return(NA)
    }
    
  }))
  
  return(tibble(chr=chromo,hub_5kb=tmp_hub_set$node,top_hub=compound_hubs))
  
}

#-----------------------------------------
union_file<-"./data/DAGGER_tbl/GM12878_union_dagger_tbl.Rda"
union_pval_file<-"./data/pval_tbl/CAGE_union_GM12878_pval_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"

dagger_hub_tbl<-input_data_fn(union_file)
cl_tbl<-input_data_fn(union_pval_file)

dagger_hub_tbl<-dagger_hub_tbl %>% 
  left_join(.,cl_tbl,by=c("node"="cl","chr"="chr","res"="res","emp.pval"="emp.pval"))


tmp_res<-"5kb"
chr_set<-unique(dagger_hub_tbl$chr)
compound_hub_l<-vector("list",length(chr_set))
names(compound_hub_l)<-chr_set


for(chromo in unique(dagger_hub_tbl$chr)){
  
  compound_hub_l[[chromo]]<-get_tophub_fn(dagger_hub_tbl,spec_res_file,chromo,tmp_res)
  
}
compound_hub_tbl<-do.call(bind_rows,compound_hub_l)

compound_hub_tbl<-compound_hub_tbl %>% 
  mutate(top_res=map_chr(top_hub,function(x){
    strsplit(x,split="_")[[1]][1]
  })) 

compound_hub_tbl %>% 
  group_by(top_res) %>% 
  summarise(n=n()) %>% 
  mutate(top_res=fct_relevel(top_res,names(res_num))) %>% 
  ggplot(.,aes(top_res,n))+geom_bar(stat="identity")

save(compound_hub_tbl,file="./data/DAGGER_tbl/GM12878_mres_coumpund_hub.Rda")
compound_hub_tbl %>% 
  filter(top_res %in% c("1Mb","500kb"))
dagger_hub_tbl %>% 
  filter(res %in% c("5kb")) %>% 
  distinct(chr,node)
