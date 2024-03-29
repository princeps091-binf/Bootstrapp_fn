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
#  filter(uni.pval<=0.5) %>% 
  ggplot(.,aes(-log10(tss.pval),-log10(uni.pval)))+
  geom_point(alpha=0.2)+
  geom_smooth()+
  # geom_density_2d_filled()+
  facet_wrap(res~.)
ggsave("~/Documents/multires_bhicect/weeklies/weekly47/img/tss_vs_uni_pval.png")

all_pval_tbl %>% filter(enh.pval<=0.0001 & tss.pval<=0.0001) %>% 
  group_by(res) %>% summarise(n=n())

all_pval_tbl %>% filter(enh.pval<=0.0001 & tss.pval<=0.0001) %>% 
  arrange(desc(uni.pval))

all_pval_tbl %>% filter(enh.pval<=0.0001 & tss.pval<=0.0001) %>% 
  filter(res=="5kb")  

enh_file<-"./data/DAGGER_tbl/HMEC_enh_dagger_tbl.Rda"
enh_dagger_tbl<-get(load(enh_file))
tmp_obj<-names(mget(load(enh_file)))
rm(list=tmp_obj)
rm(tmp_obj)


cl_tss_tbl %>% left_join(.,all_pval_tbl) %>% 
  mutate(tss.n=ifelse(is.na(tss.n),0,tss.n),enh.n=ifelse(is.na(enh.n),0,enh.n)) %>% 
  mutate(res=fct_relevel(res,res_set)) %>% 
  #  filter(uni.pval<=0.5) %>% 
  ggplot(.,aes(enh.n,tss.n))+ #scale_x_log10()+ scale_y_log10()+
  geom_point(alpha=0.2)+
#  geom_smooth()+ 
  # geom_density_2d_filled()+
  facet_wrap(res~.,scales="free")
ggsave("~/Documents/multires_bhicect/weeklies/weekly48/img/tss_n_vs_enh_n_tsshub.png")

#-----------------------------------------
## Examine distributional property of highlighted hubs

cl_union_tbl %>% 
  group_by(res) %>% summarise(n=n()) %>% 
  mutate(res=fct_relevel(res,names(res_num))) %>% 
  ggplot(.,aes(res,n))+geom_bar(stat="identity")

cl_union_tbl %>% mutate(foot=pmap_dbl(list(res,bins),function(res,bins){
  res_num[res]*length(bins)
}),
span= map_dbl(bins,function(x){
  diff(range(as.numeric(x)))
}), res=fct_relevel(res,names(res_num))) %>% 
  ggplot(.,aes(foot,span))+
  geom_point(alpha=0.05)+
#  geom_density_2d()+
  facet_wrap(res~.,scales="free")
