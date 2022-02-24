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

parent_hub_tbl<-compound_hub_tbl %>%
  filter(!(is.na(parent.hub))) %>% 
  distinct(chr,parent.hub)

parent_hub_tbl<-parent_hub_tbl %>% 
  left_join(.,cl_tbl,by=c("parent.hub"="cl","chr"="chr"))

plan(multisession, workers=3)
parent_hub_tbl<-parent_hub_tbl %>% 
#  dplyr::slice(1:30) %>% 
  mutate(hub.5kb.io=future_map_dbl(GRange,function(x){
    tmp_vec<-tss_Grange@elementMetadata$hub.5kb.io[unique(queryHits(findOverlaps(tss_Grange,x)))]
    return(sum(tmp_vec=="in")/length(tmp_vec))
  })) %>% arrange(hub.5kb.io)

parent_hub_tbl %>% 
  filter(hub.5kb.io>0.5) %>% 
  group_by(res) %>% 
  summarise(n=n())

parent_hub_tbl %>% 
  filter(hub.5kb.io>0.5 & res == "500kb") %>% 
  group_by(res) %>% 
  summarise(n=n())
#--------------------------------------
#5kb stubs
length(unique(queryHits(findOverlaps(tss_Grange,GenomicRanges::reduce(compound_hub_tbl %>% filter(is.na(parent.hub))%>% 
  left_join(.,cl_tbl,by=c("hub.5kb"="cl","chr"="chr")) %>% 
  dplyr::select(GRange) %>% unlist %>% GRangesList %>% unlist)))))

compound_hub_tbl %>% 
  distinct(chr,hub.5kb)

compound_hub_tbl %>% 
  filter(!(is.na(parent.hub))) %>% 
  group_by(chr,hub.5kb) %>% 
  summarise(n=n())

compound_hub_tbl %>% 
  filter(is.na(parent.hub)) %>% 
  distinct(chr,hub.5kb)

  