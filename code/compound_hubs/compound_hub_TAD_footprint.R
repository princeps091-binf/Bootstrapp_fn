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
union_hub_file<-"./data/pval_tbl/CAGE_union_H1_TAD_pval_tbl.Rda"
compound_hub_5kb_file<-"./data/candidate_compound_hub/H1_5kb_tss_compound_hub.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/H1/Dekker/spec_res/"

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


hmec_tad_tbl<-tbl_in_fn(union_hub_file)
hmec_tad_tbl<-do.call(bind_rows,hmec_tad_tbl %>% 
                        group_by(chr) %>% 
                        mutate(FDR=p.adjust(emp.pval,method="fdr")) %>% 
                        filter(FDR<=0.01) %>% 
                        mutate(coord=pmap(list(ID,GRange),function(ID,Grange){
                          inter_tbl<-tibble(as.data.frame(Grange)) %>% mutate(tad=ID)
                          return(inter_tbl)
                          
                        })) %>% 
                        ungroup %>% 
                        dplyr::select(coord) %>% as.list)
full_set_tbl<-hmec_tad_tbl %>% 
  dplyr::select(-tad) %>% 
  mutate(res="TAD",set='TAD') %>% 
  bind_rows(.,hmec_fdr_tbl %>% 
              mutate(set="hub"))

gg_foot<-full_set_tbl%>%
  mutate(seqnames=fct_relevel(seqnames,paste0('chr',1:22)),
         res=fct_relevel(res,c(res_set,"TAD")))%>%
  ggplot(.,aes(xmin=start,xmax=end,ymin=0,ymax=1,fill=res))+
  geom_rect()+
  facet_grid(seqnames~set)+
  scale_fill_brewer(palette = "Dark2")

gg_foot<-gg_foot +theme(axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank())+
  scale_x_continuous(labels= label_number(scale = 1/1e6,suffix="Mb"))

gg_foot
ggsave("~/Documents/multires_bhicect/weeklies/weekly55/img/H1_TAD_cl_footprint.png",height=25,width=30,units = "cm")
