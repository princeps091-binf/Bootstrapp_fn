library(tidyverse)
library(GenomicRanges)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
enh_file<-"./data/pval_tbl/CAGE_enh_HMEC_pval_tbl.Rda"
tss_file<-"./data/pval_tbl/CAGE_tss_HMEC_pval_tbl.Rda"
union_file<-"./data/pval_tbl/CAGE_union_HMEC_pval_tbl.Rda"

cl_enh_tbl<-get(load(enh_file))
tmp_obj<-names(mget(load(enh_file)))
rm(list=tmp_obj)
rm(tmp_obj)
cl_enh_tbl<-cl_enh_tbl %>% dplyr::rename(enh.pval=emp.pval,enh.n=feature_n)

cl_tss_tbl<-get(load(tss_file))
tmp_obj<-names(mget(load(tss_file)))
rm(list=tmp_obj)
rm(tmp_obj)
cl_tss_tbl<-cl_tss_tbl %>% dplyr::rename(tss.pval=emp.pval,tss.n=feature_n)

cl_union_tbl<-get(load(union_file))
tmp_obj<-names(mget(load(union_file)))
rm(list=tmp_obj)
rm(tmp_obj)
cl_union_tbl<-cl_union_tbl %>% dplyr::rename(uni.pval=emp.pval,uni.n=feature_n)

all_pval_tbl<-cl_union_tbl %>% select(-GRange) %>% full_join(.,cl_tss_tbl%>% select(-GRange)) %>% full_join(.,cl_enh_tbl%>% select(-GRange))

all_pval_tbl %>% 
  mutate(res=fct_relevel(res,res_set)) %>% 
  filter(uni.pval<=0.5) %>% 
  ggplot(.,aes(-log10(tss.pval),-log10(uni.pval)))+
  geom_point(alpha=0.1)+
  geom_smooth()+
  # geom_density_2d_filled()+
  facet_wrap(res~.)

all_pval_tbl %>% filter(enh.pval<=0.0001 & tss.pval<=0.0001) %>% 
  group_by(res) %>% summarise(n=n())

all_pval_tbl %>% filter(enh.pval<=0.0001 & tss.pval<=0.0001) %>% 
  arrange(desc(uni.pval))

all_pval_tbl %>% filter(enh.pval<=0.0001 & tss.pval<=0.0001) %>% 
  filter(res=="5kb")  
