library(tidyverse)
library(furrr)
library(data.tree)
library(igraph)
library(RColorBrewer)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
input_data_fn<-function(tmp_file){
  out_tbl<-get(base::load(tmp_file))
  tmp_obj<-names(mget(base::load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl) 
}

#-----------------------------------------
union_file<-"./data/pval_tbl/CAGE_union_HMEC_pval_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
compound_hub_5kb_file<-"./data/DAGGER_tbl/trans_res/HMEC_union_trans_res_dagger_tbl.Rda"

dagger_hub_tbl<-input_data_fn(union_file)
compound_hub_5kb_tbl<-input_data_fn(compound_hub_5kb_file)

chromo<-"chr19"
#tmp_res<-"5kb"
tmp_hub_set<-dagger_hub_tbl %>% filter(chr==chromo)
tmp_compound_hub_tbl<-compound_hub_5kb_tbl %>% filter(chr==chromo )


base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
chr_bpt<-FromListSimple(chr_spec_res$part_tree)
node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
node_ancestor<-lapply(node_ancestor,'[',-1)

#recover all CAGE containing leaves
chr_hub_set<-unique(c(tmp_hub_set$cl,unique(unlist(node_ancestor[tmp_hub_set$cl])))) 
#rebuild corresponding tree
Prune(chr_bpt, function(x) x$name %in% chr_hub_set)
g<-graph_from_data_frame(ToDataFrameNetwork(chr_bpt))
chr_hubs<-dagger_hub_tbl %>% filter(chr==chromo ) %>% distinct(cl) %>% unlist

# Plot BPTs with color code of choice
cl_res_vec<-tmp_hub_set$res
names(cl_res_vec)<-tmp_hub_set$cl
cl_size<- rep(0.5,length(V(g)$name))
names(cl_size)<-V(g)$name
in_set<-tmp_hub_set %>% filter(emp.pval<= 1e-4) %>% 
  dplyr::select(cl) %>% unlist
cl_size[in_set]<-2

res_col_pal<-brewer.pal(6,"Dark2")
names(res_col_pal)<-res_set


V(g)$color<-res_col_pal[cl_res_vec[V(g)$name]]
V(g)$size<-cl_size

plot(g,layout=layout_as_tree(g),vertex.frame.color=NA,vertex.size=2,vertex.label=NA,edge.arrow.size=0)
legend('topright',legend=res_set,col=res_col_pal,pch=16)

plot(g,layout=layout_as_tree(g),vertex.frame.color=NA,vertex.label=NA,edge.arrow.size=0)
legend('topright',legend=res_set,col=res_col_pal,pch=16)

png(filename = "~/Documents/multires_bhicect/weeklies/group_meeting/group_meeting_04_2022/img/HMEC_chr19_all_cage_clusters_bpt.png",width=20,height = 20,units = "cm",res=500)
par(mar=c(0,0,0,0))
plot(g,layout=layout_as_tree(g),vertex.frame.color=NA,vertex.size=2,vertex.label=NA,edge.arrow.size=0)
legend('topright',legend=res_set,col=res_col_pal,pch=16)
dev.off()

png(filename = "~/Documents/multires_bhicect/weeklies/group_meeting/group_meeting_04_2022/img/HMEC_chr19_all_cage_clusters_pval_bpt.png",width=20,height = 20,units = "cm",res=500)
par(mar=c(0,0,0,0))
plot(g,layout=layout_as_tree(g),vertex.frame.color=NA,vertex.label=NA,edge.arrow.size=0)
legend('topright',legend=res_set,col=res_col_pal,pch=16)
dev.off()


compound_set<-unique(tmp_compound_hub_tbl$node)

col_pal<-brewer.pal(3,"Set1")

tmp_vec<-cl_res_vec
tmp_vec[!(names(tmp_vec) %in% compound_set)]<-NA
tmp_col_vec<-res_col_pal[tmp_vec[V(g)$name]]
tmp_col_vec[is.na(tmp_col_vec)]<-"grey70"

V(g)$color<-tmp_col_vec
V(g)$size<-ifelse(V(g)$name %in% compound_set,2,0.5)

plot(g,layout=layout_as_tree(g),vertex.frame.color=NA,vertex.label=NA,edge.arrow.size=0)
png(filename = "~/Documents/multires_bhicect/weeklies/group_meeting/group_meeting_04_2022/img/HMEC_chr19_all_DAGGER_hubs_bpt.png",width=20,height = 20,units = "cm",res=500)
par(mar=c(0,0,0,0))
plot(g,layout=layout_as_tree(g),vertex.frame.color=NA,vertex.label=NA,edge.arrow.size=0)
legend('topright',legend=res_set,col=res_col_pal,pch=16)
dev.off()
#---------------------
dagger_hub_tbl %>% 
  mutate(res=fct_relevel(res,res_set)) %>% 
  ggplot(.,aes(emp.pval))+
  geom_histogram()+
  facet_wrap(res~.,scales="free")
ggsave("~/Documents/multires_bhicect/weeklies/group_meeting/group_meeting_04_2022/img/HMEC_pval_hist.png")
