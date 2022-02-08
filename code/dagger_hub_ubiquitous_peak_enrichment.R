library(GenomicRanges)
library(tidyverse)
library(furrr)
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

build_Grange<-function(tmp_tbl){
  
  cage_chr_Grange<-GRanges(seqnames=tmp_tbl$chr,
                           ranges = IRanges(start=as.numeric(tmp_tbl$start),
                                            end=as.numeric(tmp_tbl$end)
                           ))
  return(cage_chr_Grange)
  
}
#-----------------------------------------
union_file<-"./data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"
peak_summary_file<-"~/Documents/multires_bhicect/Fantom5_CAGE_peak_content_analysis/data/CAGE_tss_cell_line_summary_tbl.Rda"
HMEC_peak_file<-"./data/CAGE_tss_coord_HMEC_tbl.Rda"
union_pval_file<-"./data/pval_tbl/CAGE_union_HMEC_pval_tbl.Rda"


dagger_hub_tbl<-input_data_fn(union_file)
cl_tbl<-input_data_fn(union_pval_file)

peak_summary_tbl<-input_data_fn(peak_summary_file)
HMEC_peak_summary_tbl<-input_data_fn(HMEC_peak_file)

peak_Grange<-build_Grange(peak_summary_tbl %>% filter(peak.ID %in% HMEC_peak_summary_tbl$Id))
mcols(peak_Grange)<-peak_summary_tbl %>% filter(peak.ID %in% HMEC_peak_summary_tbl$Id) %>% select(peak.ID,n,mad,med)

plan(multisession, workers = 3)

cl_tbl<-cl_tbl %>% 
#  slice(1:5) %>% 
  mutate(enh.sum=future_map(GRange,function(x){
  as_tibble(peak_Grange@elementMetadata) %>% slice(unique(queryHits(findOverlaps(peak_Grange,x))))
}))

hub_enh_tbl<-dagger_hub_tbl %>% filter(res=="500kb") %>% 
  left_join(.,cl_tbl,by=c("node"="cl","chr"="chr","res"="res","emp.pval"="emp.pval")) %>% 
  select(enh.sum) %>% 
  as.list() %>% do.call(bind_rows,.) %>% 
  distinct 



hub_enh_tbl %>% mutate(set="hub") %>% 
  bind_rows(.,peak_summary_tbl %>% filter(peak.ID %in% HMEC_peak_summary_tbl$Id) %>% 
              filter(!(peak.ID %in% hub_enh_tbl$peak.ID)) %>% 
              mutate(set="out")) %>% 
  ggplot(.,aes(n,color=set))+geom_density()+
  geom_density(data=peak_summary_tbl %>% filter(peak.ID %in% HMEC_peak_summary_tbl$Id) %>% mutate(set="all"),aes(n))
ggsave("~/Documents/multires_bhicect/weeklies/weekly50/img/HMEC_500kb_hub_vs_out_tss_ubiquity.png")
