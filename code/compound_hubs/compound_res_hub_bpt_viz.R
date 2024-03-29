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

#-----------------------------------------
union_file<-"./data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"
compound_cl_file<-"./data/HMEC_5kb_hub_ancestry_5kb_peak_compound_hub.Rda"
hub_peak_content_tbl_file<-"./data/HMEC_5kb_hub_ancestry_CAGE_peak_content.Rda"
compound_hub_tbl_file<-"./data/HMEC_5kb_hub_ancestry.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"

dagger_hub_tbl<-input_data_fn(union_file)
compound_cl<-input_data_fn(compound_cl_file)
hub_peak_content_tbl<-input_data_fn(hub_peak_content_tbl_file)
compound_hub_tbl<-input_data_fn(compound_hub_tbl_file)

chromo<-"chr14"
tmp_res<-"5kb"
tmp_hub_set<-dagger_hub_tbl %>% filter(chr==chromo)

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

png(filename = "~/Documents/multires_bhicect/weeklies/weekly52/img/chr14_BPT_dagger_hub.png")
plot(g,layout=layout_as_tree(g),vertex.size=2,vertex.label=NA,edge.arrow.size=0)
dev.off()
chr_compound_hubs<-compound_cl %>% filter(chr==chromo & !(grepl("5kb_",hub))) %>% dplyr::select(hub) %>% unlist

V(g)$color<-ifelse(V(g)$name %in% chr_compound_hubs,"green",ifelse(V(g)$name %in% chr_hubs,"red","grey50"))
plot(g,layout=layout_as_tree(g),vertex.size=3,vertex.label=NA,edge.arrow.size=0)
png(filename = "~/Documents/multires_bhicect/weeklies/weekly52/img/chr14_BPT_compound_hub.png")
plot(g,layout=layout_as_tree(g),vertex.size=2,vertex.label=NA,edge.arrow.size=0)
dev.off()
