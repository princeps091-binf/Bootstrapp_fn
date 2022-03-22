library(tidyverse)
library(scales)
library(GenomicRanges)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------------------------------------------
tbl_in_fn<-function(tmp_file){
  tmp_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  
  return(tmp_tbl)
}
#-------------------------------------------------------------------------------------------------------

union_hub_file<-"./data/pval_tbl/CAGE_union_H1_TAD_pval_tbl.Rda"

hmec_dagger_01_tbl<-tbl_in_fn(union_hub_file)
hmec_fdr_tbl<-do.call(bind_rows,hmec_dagger_01_tbl %>% 
  group_by(chr) %>% 
  mutate(FDR=p.adjust(emp.pval,method="fdr")) %>% 
  filter(FDR<=0.01) %>% 
  mutate(coord=pmap(list(ID,GRange),function(ID,Grange){
    inter_tbl<-tibble(as.data.frame(Grange)) %>% mutate(tad=ID)
    return(inter_tbl)
    
  })) %>% 
    ungroup %>% 
  dplyr::select(coord) %>% as.list)

gg_foot<-hmec_fdr_tbl%>%
  mutate(seqnames=fct_relevel(seqnames,paste0('chr',1:22)))%>%
  ggplot(.,aes(xmin=start,xmax=end,ymin=0,ymax=1))+
  geom_rect()+
  facet_grid(seqnames~.)+
  scale_fill_brewer(palette = "Dark2")

gg_foot<-gg_foot +theme(axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        legend.position="none")+
  scale_x_continuous(labels= label_number(scale = 1/1e6,suffix="Mb"))

gg_foot
ggsave(filename = "~/Documents/multires_bhicect/weeklies/weekly53/img/H1_TAD_foot.png",height = 25,width=15,units = "cm",gg_foot)
