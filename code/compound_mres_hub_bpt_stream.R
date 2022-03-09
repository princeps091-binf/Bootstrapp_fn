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
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl) 
}

#-----------------------------------------
union_file<-"./data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
compound_hub_5kb_file<-"./data/candidate_compound_hub/HMEC_5kb_tss_compound_hub.Rda"

dagger_hub_tbl<-input_data_fn(union_file)
compound_hub_5kb_tbl<-input_data_fn(compound_hub_5kb_file)

chromo<-"chr19"
tmp_res<-"5kb"
tmp_hub_set<-dagger_hub_tbl %>% filter(chr==chromo & res==tmp_res)
tmp_compound_hub_tbl<-compound_hub_5kb_tbl %>% filter(chr==chromo )


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

top_compound_hubs<-tmp_compound_hub_tbl %>% 
  filter(parent.hub %in% tmp_compound_hub_tbl$parent.hub[which(!(tmp_compound_hub_tbl$parent.hub %in% unique(unlist(tmp_compound_hub_tbl$ch.hub))))]) %>% 
  dplyr::select(parent.hub) %>% unlist

compound_set<-unique(c(tmp_compound_hub_tbl$parent.hub,unlist(tmp_compound_hub_tbl$ch.hub)))

V(g)$color<-ifelse(V(g)$name %in% top_compound_hubs,"green",
                   ifelse(V(g)$name %in% compound_set,"orange",
                          ifelse(V(g)$name %in% chr_hubs, "red","grey50")))
plot(g,layout=layout_as_tree(g),vertex.size=2,vertex.label=NA,edge.arrow.size=0)


