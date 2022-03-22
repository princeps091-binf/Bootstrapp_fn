## footprint stream
library(tidyverse)
library(scales)
library(GenomicRanges)
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
cl_reduce_coord_fn<-function(hmec_dagger_01_tbl,tmp_res,res_num){
  
  tmp_bin_tbl<-hmec_dagger_01_tbl %>% filter(res==tmp_res) %>%unnest(cols = "bins")
  tmp_bin_tbl<-tmp_bin_tbl%>%mutate(end = as.numeric(bins) + res_num[res]-1)
  tmp_bin_tbl<-tmp_bin_tbl%>%distinct(chr,bins,end)
  inter_cl_Grange<-   GRanges(seqnames=tmp_bin_tbl$chr,
                              ranges = IRanges(start=as.numeric(tmp_bin_tbl$bins),
                                               end=tmp_bin_tbl$end,
                                               names=paste("cl_inter",1:nrow(tmp_bin_tbl),sep='_')
                              ))
  inter_cl_Grange<-GenomicRanges::reduce(inter_cl_Grange)
  inter_tbl<-tibble(as.data.frame(inter_cl_Grange))
  return(inter_tbl)
  
}

#-------------------------------------------------------------------------------------------------------
compound_hub_5kb_file<-"./data/candidate_compound_hub/GM12878_5kb_tss_compound_hub.Rda"
union_cl_file<-"./data/pval_tbl/CAGE_union_GM12878_pval_tbl.Rda"

compound_hub_5kb_tbl<-tbl_in_fn(compound_hub_5kb_file)
cl_union_tbl<-tbl_in_fn(union_cl_file)

hmec_dagger_01_tbl<-compound_hub_5kb_tbl %>% left_join(.,cl_union_tbl %>% dplyr::rename(parent.hub=cl)%>% dplyr::select(chr,parent.hub,res,bins))


hmec_fdr_tbl<-do.call(bind_rows,lapply(unique(hmec_dagger_01_tbl$res),function(f){
  return(cl_reduce_coord_fn(hmec_dagger_01_tbl,f,res_num)%>%mutate(res=f))
}))


gg_foot<-hmec_fdr_tbl%>%
  mutate(seqnames=fct_relevel(seqnames,paste0('chr',1:22)),res=fct_relevel(res,c(res_set[res_set %in% .$res])))%>%
  ggplot(.,aes(xmin=start,xmax=end,ymin=0,ymax=1,fill=as.factor(res)))+
  geom_rect()+
  facet_grid(seqnames~.)+
  scale_fill_brewer(palette = "Set1")

gg_foot<-gg_foot +theme(axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank())+
  scale_x_continuous(labels= label_number(scale = 1/1e6,suffix="Mb"))

gg_foot
