library(data.tree)
library(tidyverse)
library(furrr)

options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------
feature_pval_file<-"./data/pval_tbl/CAGE_union_GM12878_pval_tbl.Rda"
BHiCect_res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"

feature_pval_tbl<-get(load(feature_pval_file))
tmp_obj<-names(mget(load(feature_pval_file)))
rm(list=tmp_obj)
rm(tmp_obj)
chr_set<-unique(feature_pval_tbl$chr)

pval_inheritance_l<-lapply(chr_set,function(chromo){
  chr_pval_tbl<-feature_pval_tbl %>% filter(chr==chromo)
  cat(chromo, ": select clusters with two feature-containing bins \n")
  load(paste0(BHiCect_res_file,chromo,"_spec_res.Rda"))
  cat(chromo, ": Build tree \n")
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  #build the cage-containing sub-tree
  cage_node<-unlist(chr_pval_tbl%>%dplyr::select(cl))
  cage_set<-unique(c(cage_node,unique(unlist(node_ancestor[cage_node])))) 
  Prune(chr_bpt, function(x) x$name %in% cage_set)
  p_node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  
  #rebuild corresponding tree according to DAGGER
  node_parent<-lapply(p_node_ancestor,'[',2)
  cl_pval<-chr_pval_tbl$emp.pval
  names(cl_pval)<-chr_pval_tbl$cl
  
  return(tibble(child.pval=cl_pval[names(node_parent)],parent.pval=cl_pval[unlist(node_parent)],parent=unlist(node_parent),child=names(node_parent)) %>% 
           mutate(parent.res=str_split_fixed(parent,"_",2)[,1],child.res=str_split_fixed(child,"_",2)[,1])) 
})
do.call(bind_rows,pval_inheritance_l) %>% 
  filter(parent.res==child.res) %>%
  ggplot(.,aes(log10(child.pval),log10(parent.pval)))+
  geom_point(alpha=0.01)
