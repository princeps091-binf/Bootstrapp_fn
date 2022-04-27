library(tidyverse)
library(furrr)
library(data.tree)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
get_obj_in_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#-----------------------------------------
# explicit evaluation of ancestry
spec_res_file<-"~/Documents/multires_bhicect/data/H1/Dekker/spec_res/"
hub_file<-"./data/DAGGER_tbl/trans_res/H1_union_trans_res_dagger_tbl.Rda"
out_file<-"./data/DAGGER_tbl/trans_res/H1_union_top_trans_res_dagger_tbl.Rda"

hub_tbl<-get_obj_in_fn(hub_file)
top_cl_tbl<-do.call(bind_rows,map(unique(hub_tbl$chr),function(chromo){
  message(chromo)
  base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  tmp_tbl<-hub_tbl %>% 
    filter(chr==chromo)
  cl_set<-tmp_tbl$node[unlist(lapply(tmp_tbl$node,function(x){
    !(any(tmp_tbl$node %in% node_ancestor[[x]]))
  }))]  
  return(tibble(chr=chromo,node=cl_set))
}))

top_cl_tbl %>% 
  left_join(.,hub_tbl) %>% 
  group_by(res) %>% 
  summarise(n=n()) %>% 
  arrange(n)
save(top_cl_tbl,file=out_file)
