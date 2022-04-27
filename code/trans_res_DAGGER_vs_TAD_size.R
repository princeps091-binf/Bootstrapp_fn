# Top-hub generic topo stat
library(tidyverse)
library(GenomicRanges)
library(scales)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------------------------------------------
tbl_in_fn<-function(tmp_file){
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  
  return(out_tbl)
}

Build_coord_fn<-function(top_compound_hub_5kb_tbl,spec_res_file){
  coord_tbl<-do.call(bind_rows,map(unique(top_compound_hub_5kb_tbl$chr),function(chromo){
    message(chromo)
    base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
    tmp_tbl<-top_compound_hub_5kb_tbl %>% 
      filter(chr==chromo) %>% 
      mutate(bins=chr_spec_res$cl_member[node]) %>% 
      mutate(bins=map(bins,as.numeric)) 
    
  }))
}

Build_GRange_fn<-function(chromo,res,bins,res_num){
  inter_cl_Grange<-   GRanges(seqnames=chromo,
                              ranges = IRanges(start=bins,
                                               end=bins + res_num[res]-1
                              ))
  inter_cl_Grange<-GenomicRanges::reduce(inter_cl_Grange)
  return(inter_cl_Grange)
  
}

topo_summary_fn<-function(top_compound_hub_5kb_tbl,spec_res_file,cage_GRange){
  message("Build GRange")
  top_compound_hub_5kb_tbl<-Build_coord_fn(top_compound_hub_5kb_tbl,spec_res_file) %>% 
    mutate(GRange=pmap(list(chr,res,bins),function(chromo,res,bins){
      Build_GRange_fn(chromo,res,bins,res_num)
    }))

  message("Compute topo stats")
  
  tmp_tbl<-top_compound_hub_5kb_tbl %>% 
    mutate(foot=map_int(GRange,function(x)sum(width(x)))) %>%
    mutate(res=fct_relevel(res,names(res_num))) %>% 
    dplyr::select(chr,res,node,foot)
  return(tmp_tbl)
  
}
#-------------------------------------------------------------------------------------------------------
hub_files<-list(H1="./data/DAGGER_tbl/trans_res/H1_union_top_trans_res_dagger_tbl.Rda",
                GM12878="./data/DAGGER_tbl/trans_res/GM12878_union_top_trans_res_dagger_tbl.Rda",
                HMEC="./data/DAGGER_tbl/trans_res/HMEC_union_top_trans_res_dagger_tbl.Rda")

spec_res_files<-list(H1="~/Documents/multires_bhicect/data/H1/Dekker/spec_res/",
                     GM12878="~/Documents/multires_bhicect/data/GM12878/spec_res/",
                     HMEC="~/Documents/multires_bhicect/data/HMEC/spec_res/")

TAD_files<-list(H1="./data/pval_tbl/CAGE_union_H1_TAD_pval_tbl.Rda",
                GM12878="./data/pval_tbl/CAGE_union_GM12878_TAD_pval_tbl.Rda",
                HMEC="./data/pval_tbl/CAGE_union_HMEC_TAD_pval_tbl.Rda")
#-------------------------------------------------------------------------------------------------------
## Subset the mutually exclusive collection of compound hubs
topo_summary_tbl<-do.call(bind_rows,lapply(1:length(hub_files),function(cell_line){

  TAD_tbl<-tbl_in_fn(TAD_files[[cell_line]]) %>% 
    mutate(foot=map_int(GRange,function(x){
      sum(width(x))
    }),
    cell.line=names(spec_res_files)[cell_line],
    res="TAD",set="TAD") %>% 
    dplyr::rename(node=ID) %>% 
    dplyr::select(chr,res,node,foot,set,cell.line)
  
  hub_tbl<-tbl_in_fn(hub_files[[cell_line]]) %>% 
    mutate(res=str_split_fixed(node,"_",2)[,1])
  
  tmp_tbl<-topo_summary_fn(hub_tbl,spec_res_files[[cell_line]],cage_GRange) %>% 
    mutate(cell.line=names(spec_res_files)[cell_line],set="hubs")
  tmp_tbl<-tmp_tbl %>% 
    bind_rows(.,TAD_tbl)
  return(tmp_tbl)
}))

topo_summary_tbl %>% 
  filter(set!="TAD") %>% 
  bind_rows(.,
  do.call(bind_rows,lapply(unique(topo_summary_tbl$res),function(r){
  topo_summary_tbl %>% filter(set=="TAD") %>% 
    dplyr::select(-res) %>% 
    mutate(res=r)
  
  }))) %>% 
  mutate(res=fct_relevel(res,res_set)) %>% 
  filter(res!="TAD") %>% 
  ggplot(.,aes(foot,color=set))+
  geom_density()+
  scale_x_log10()+
  facet_grid(res~cell.line,scales="free")+
  scale_color_brewer(palette="Set1")
