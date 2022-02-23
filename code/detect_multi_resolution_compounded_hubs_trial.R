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

get_top_cl_fn<-function(cl_tbl,res_file){
  chr_set<-unique(cl_tbl$chr)
  res_l<-vector("list",length(chr_set))
  names(res_l)<-chr_set
  for(chromo in unique(cl_tbl$chr)){
    message(chromo)
    chr_cl_tbl<-cl_tbl%>%filter(chr==chromo)
    base::load(paste0(res_file,chromo,"_spec_res.Rda"))
    chr_bpt<-FromListSimple(chr_spec_res$part_tree)
    node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
    node_ancestor<-lapply(node_ancestor,'[',-1)
    res_l[[chromo]]<-tibble(chr=chromo,top.cl=chr_cl_tbl$node[which(unlist(lapply(chr_cl_tbl$node,function(x){
      !(any(chr_cl_tbl$node %in% node_ancestor[[x]]))
    })))])
  }  
  return(do.call(bind_rows,res_l))
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
union_file<-"./data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"
union_pval_file<-"./data/pval_tbl/CAGE_union_HMEC_pval_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"

dagger_hub_tbl<-input_data_fn(union_file)
cl_tbl<-input_data_fn(union_pval_file)

dagger_hub_tbl<-dagger_hub_tbl %>% 
  left_join(.,cl_tbl,by=c("node"="cl","chr"="chr","res"="res","emp.pval"="emp.pval"))

dagger_hub_tbl %>% 
  mutate(res=fct_relevel(res,names(res_num))) %>% 
#  mutate(foot=map_int(GRange,function(x)sum(width(x)))) %>%
  mutate(nbin=map_int(bins,function(x)length(x))) %>% 
  ggplot(.,aes(foot,color=res))+geom_density()+facet_wrap(res~.,scales="free")


dagger_hub_tbl %>% 
  mutate(res=fct_relevel(res,names(res_num))) %>% 
  #  mutate(foot=map_int(GRange,function(x)sum(width(x)))) %>%
  mutate(nbin=map_int(bins,function(x)length(x))) %>% 
  filter(nbin==2) %>% 
  mutate(span=unlist(pmap(list(res,bins),function(res,bins){
    diff(range(as.numeric(bins)))/res_num[res]
  }))) %>% group_by(res,span) %>% 
  summarise(n=n()) %>% 
  filter(res=="5kb") %>% 
  arrange(desc(n))


top_hub_5kb<-get_top_cl_fn(dagger_hub_tbl %>% filter(res=="5kb"),spec_res_file)

top_hub_5kb %>% 
  left_join(.,dagger_hub_tbl,by=c("top.cl"="node","chr"='chr')) %>% 
  mutate(span=unlist(pmap(list(res,bins),function(res,bins){
    diff(range(as.numeric(bins)))/res_num[res]
  }))) 
#-----------------------------------------------
# Examine BPT distribution of DAGGER hubs
chromo<-"chr19"
tmp_res<-"5kb"
tmp_hub_set<-dagger_hub_tbl %>% filter(chr==chromo & res==tmp_res)

base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
chr_bpt<-FromListSimple(chr_spec_res$part_tree)
node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
node_ancestor<-lapply(node_ancestor,'[',-1)

#recover all CAGE containing leaves
chr_hub_set<-unique(c(tmp_hub_set$node,unique(unlist(node_ancestor[tmp_hub_set$node])))) 
#rebuild corresponding tree
Prune(chr_bpt, function(x) x$name %in% chr_hub_set)
g<-graph_from_data_frame(ToDataFrameNetwork(chr_bpt))
chr_hubs<-dagger_hub_tbl %>% filter(chr==chromo ) %>% distinct(node) %>% unlist
V(g)$color<-ifelse(V(g)$name %in% chr_hubs,"red","grey50")
plot(g,layout=layout_as_tree(g),vertex.size=2,vertex.label=NA,edge.arrow.size=0)

compound_hubs<-unique(unlist(lapply(node_ancestor[tmp_hub_set$node],function(x){
  x[max(which(x %in% chr_hubs))]
})))
V(g)$color<-ifelse(V(g)$name %in% tmp_hub_set$node,"red",ifelse(V(g)$name %in% compound_hubs, "green",ifelse(V(g)$name %in% chr_hubs,"orange","grey50")))
plot(g,layout=layout_as_tree(g),vertex.size=2,vertex.label=NA,edge.arrow.size=0)
#-----------------------------------------------

tmp_res<-"5kb"
chr_set<-unique(dagger_hub_tbl$chr)
compound_hub_l<-vector("list",length(chr_set))
names(compound_hub_l)<-chr_set


for(chromo in unique(dagger_hub_tbl$chr)){

  compound_hub_l[[chromo]]<-get_tophub_fn(dagger_hub_tbl,spec_res_file,chromo,tmp_res)

}
compound_hub_tbl<-do.call(bind_rows,compound_hub_l)
compound_hub_tbl %>% 
  mutate(top_res=map_chr(top_hub,function(x){
    strsplit(x,split="_")[[1]][1]
  })) %>% 
  group_by(top_res) %>% 
  summarise(n=n()) %>% 
  mutate(top_res=fct_relevel(top_res,names(res_num))) %>% 
  ggplot(.,aes(top_res,n))+geom_bar(stat="identity")

compound_hub_tbl<-compound_hub_tbl %>% 
  mutate(top_res=map_chr(top_hub,function(x){
    strsplit(x,split="_")[[1]][1]
  })) 
save(compound_hub_tbl,file="./data/DAGGER_tbl/HMEC_mres_coumpund_hub.Rda")
compound_hub_tbl %>% 
  filter(top_res %in% c("1Mb","500kb"))
dagger_hub_tbl %>% 
  filter(res %in% c("5kb")) %>% 
  distinct(chr,node)
