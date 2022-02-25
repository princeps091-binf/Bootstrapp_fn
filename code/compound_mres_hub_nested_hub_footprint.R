#Examine extent to which hub footprint carries over
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

Build_GRange_fn<-function(bin_set,res,chr,res_num){
  
  return(GRanges(seqnames=chr,
                 ranges = IRanges(start=as.numeric(bin_set),
                                  end=as.numeric(bin_set)+res_num[res]-1
                 )))
  
}

get_tophub_fn<-function(dagger_hub_tbl,spec_res_file,chromo,tmp_res){
  message(chromo)
  tmp_hub_set<-dagger_hub_tbl %>% filter(chr==chromo & res==tmp_res)
  chr_hubs<-dagger_hub_tbl %>% filter(chr==chromo ) %>% distinct(node) %>% unlist
  
  base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  compound_hubs<-do.call(bind_rows,lapply(tmp_hub_set$node,function(x){
    if(sum(node_ancestor[[x]] %in% chr_hubs)>0){
      tmp_vec<-node_ancestor[[x]]
      return(tibble(chr=chromo,hub.5kb=x,parent.hub=tmp_vec[which(tmp_vec %in% chr_hubs)]))
    } else{
      return(tibble(chr=chromo,hub.5kb=x,parent.hub=NA))
    }
    
  }))
  
  return(compound_hubs)
  
}

get_children_hub<-function(parent_hub_tbl,dagger_mres_hub_tbl,chromo,spec_res_file){
  
  message(chromo)
  tmp_hub_set<-parent_hub_tbl %>% filter(chr==chromo)
  chr_hubs<-dagger_mres_hub_tbl %>% filter(chr==chromo ) %>% distinct(node) %>% unlist
  
  base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  compound_hubs<-do.call(bind_rows,lapply(tmp_hub_set$parent.hub,function(x){
    tmp_ch<-names(which(unlist(lapply(node_ancestor[chr_hubs],function(y) x %in% y))))
    return(tibble(chr=chromo,parent.hub=x,children.hub=tmp_ch))
  }))
  
  return(compound_hubs)
  
}
#-------------------------------
tss_Grange_file<-"./data/GRanges/CAGE_tss_HMEC_Grange.Rda"
pval_tbl_file<-"./data/pval_tbl/CAGE_union_HMEC_pval_tbl.Rda"
dagger_union_file<-"./data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"

tss_Grange<-get_obj_in_fn(tss_Grange_file)
pval_tbl<-get_obj_in_fn(pval_tbl_file)
dagger_mres_hub_tbl<-get_obj_in_fn(dagger_union_file)


hires_hub_GRange<-GenomicRanges::reduce(dagger_mres_hub_tbl %>% 
                                          filter(res=="5kb") %>% 
                                          left_join(.,pval_tbl %>% 
                                                      dplyr::select(chr,cl,GRange),by=c("chr"='chr','node'="cl")) %>% 
                                          dplyr::select(GRange) %>% unlist %>% GRangesList %>% unlist)

hub_5kb_io_vec<-rep("out",length(tss_Grange))
hub_5kb_io_vec[unique(subjectHits(findOverlaps(hires_hub_GRange,tss_Grange)))]<-"in"
mcols(tss_Grange)<-tibble(hub.5kb.io=hub_5kb_io_vec)

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

save(compound_hub_tbl,file="./data/HMEC_5kb_hub_ancestry.Rda")

parent_hub_tbl<-compound_hub_tbl %>%
  filter(!(is.na(parent.hub))) %>% 
  distinct(chr,parent.hub)

parent_hub_tbl<-parent_hub_tbl %>% 
  left_join(.,pval_tbl,by=c("parent.hub"="cl","chr"="chr"))
#---------------------------------------------------------------------
## Collect for each parent-hub the nested hub content
chr_set<-unique(parent_hub_tbl$chr)
parent_hub_content_l<-vector("list",length(chr_set))
names(parent_hub_content_l)<-chr_set

for(chromo in chr_set){
  
  parent_hub_content_l[[chromo]]<-get_children_hub(parent_hub_tbl,dagger_mres_hub_tbl,chromo,spec_res_file)
  
}

parent_hub_content_tbl<-do.call(bind_rows,parent_hub_content_l)

save(parent_hub_content_tbl,file="./data/HMEC_parent_hub_content.Rda")

parent_hub_content_tbl %>% 
  group_by(chr,parent.hub) %>% 
  summarise(n=n()) %>% 
  mutate(parent.res=map_chr(parent.hub,function(x){
    
    return(strsplit(x,split="_")[[1]][1])
  })) %>% 
  mutate(parent.res=fct_relevel(parent.res,names(res_num))) %>% 
  ggplot(.,aes(n))+geom_density()+facet_wrap(parent.res~.,scales="free")

test<-parent_hub_content_tbl  %>%
  group_by(chr,parent.hub) %>% 
  summarise(ch.hub=list(unique(children.hub))) %>% 
  ungroup() %>% 
  mutate(parent.res=map_chr(parent.hub,function(x){
    
    return(strsplit(x,split="_")[[1]][1])
  })) %>% 
  filter(parent.res=="100kb") %>% 
  dplyr::slice_sample(n=10)

test %>% 
  mutate(hub.foot=pmap_dbl(list(chr,parent.hub,ch.hub),function(chromo,parent.hub,ch.hub){
    parent_GRange<-unlist(GenomicRanges::reduce(pval_tbl %>% dplyr::select(chr,cl,GRange) %>% filter(chr==chromo & cl == parent.hub) %>% 
      dplyr::select(GRange) %>% unlist %>% GRangesList))
    child_GRange<-GenomicRanges::reduce(Reduce(append,pval_tbl %>% dplyr::select(chr,cl,res,GRange) %>% filter(chr==chromo & cl%in% ch.hub ) %>% 
                                                 dplyr::select(GRange) %>% unlist))
    return(sum(width(child_GRange))/sum(width(parent_GRange)))
  })) %>% arrange(desc(hub.foot))

full_hub_set_tbl<-parent_hub_content_tbl %>% 
  group_by(chr) %>% 
  summarise(hubs=unique(c(parent.hub,children.hub)))

full_hub_GRange_tbl<-pval_tbl %>% dplyr::select(chr,cl,res,GRange) %>% 
  inner_join(.,full_hub_set_tbl,by=c("chr"="chr","cl"="hubs"))
  
tmp_tbl<-parent_hub_content_tbl  %>%
  group_by(chr,parent.hub) %>% 
  summarise(ch.hub=list(unique(children.hub))) %>% 
  ungroup() 

plan(multisession, workers = 5)

tmp_tbl<-tmp_tbl %>%
  mutate(hub.foot=future_pmap_dbl(list(chr,parent.hub,ch.hub),function(chromo,parent.hub,ch.hub){
    parent_GRange<-unlist(GenomicRanges::reduce(full_hub_GRange_tbl %>% filter(chr==chromo & cl == parent.hub) %>% 
                                                  dplyr::select(GRange) %>% unlist %>% GRangesList))
    child_GRange<-GenomicRanges::reduce(Reduce(append,full_hub_GRange_tbl %>% filter(chr==chromo & cl%in% ch.hub ) %>% 
                                                 dplyr::select(GRange) %>% unlist))
    return(sum(width(child_GRange))/sum(width(parent_GRange)))
  }))

tmp_tbl %>% 
  mutate(parent.res=map_chr(parent.hub,function(x){
    
    return(strsplit(x,split="_")[[1]][1])
  })) %>%
  mutate(parent.res=fct_relevel(parent.res,names(res_num))) %>% 
  ggplot(.,aes(hub.foot))+
  geom_density()+
  facet_wrap(parent.res~.,scales="free")
