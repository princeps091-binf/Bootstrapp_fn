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
spec_res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"

compound_hub_5kb_tbl<-tbl_in_fn(compound_hub_5kb_file)


top_compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(compound_hub_5kb_tbl$chr),function(chromo){
  tmp_tbl<-compound_hub_5kb_tbl %>% 
    filter(chr==chromo)
  return(tmp_tbl %>% 
           filter(parent.hub %in% tmp_tbl$parent.hub[which(!(tmp_tbl$parent.hub %in% unique(unlist(tmp_tbl$ch.hub))))])
  )
  
}))


top_compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(top_compound_hub_5kb_tbl$chr),function(chromo){
  message(chromo)
  base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  tmp_tbl<-top_compound_hub_5kb_tbl %>% 
    filter(chr==chromo) %>% 
    mutate(bins=chr_spec_res$cl_member[parent.hub]) %>% 
    mutate(bins=map(bins,as.numeric)) 

}))


hmec_fdr_tbl<-do.call(bind_rows,lapply(unique(top_compound_hub_5kb_tbl$res),function(f){
  return(cl_reduce_coord_fn(top_compound_hub_5kb_tbl,f,res_num)%>%mutate(res=f))
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
