library(data.tree)
library(tidyverse)
library(furrr)
library(GenomicRanges)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------
get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(base::load(tmp_file))
  tmp_obj<-names(mget(base::load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

#--------------------------
dagger_tbl_file<-"./data/DAGGER_tbl/trans_res/FDR_sensitivity/GM12878_union_trans_res_mFDR_dagger_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"
CAGE_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_union_GM12878_Grange.Rda"

CAGE_GRange<-get_tbl_in_fn(CAGE_GRange_file)

dagger_mres_tbl<-get_tbl_in_fn(dagger_tbl_file)

dagger_mres_tbl %>% 
  group_by(FDR) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(FDR,n))+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept = 0.01)

# Compute the footprint for different FDR threshold
## Extract the bin composition for every distinct candidate hubs found overall
candidate_hub_tbl<-dagger_mres_tbl %>% 
  dplyr::select(chr,node) %>% 
  distinct()

plan(multisession,workers=4)
candidate_hub_tbl<-candidate_hub_tbl %>% 
  mutate(bins=future_pmap(list(chr,node),function(chromo,cl){
    chr_spec_res<-get_tbl_in_fn(paste0(spec_res_file,chromo,"_spec_res.Rda"))
    return(as.numeric(chr_spec_res$cl_member[[cl]]))
  }))
plan(sequential)

candidate_hub_tbl<-candidate_hub_tbl %>% 
  mutate(res=str_split_fixed(node,"_",2)[,1])

fdr_cand_Grange_l<-lapply(unique(dagger_mres_tbl$FDR),function(i){
  message(i)
  tmp_hub_tbl<-dagger_mres_tbl %>% 
    filter(FDR <= i)
  
  tmp_hub_tbl<-tmp_hub_tbl %>% 
    left_join(.,candidate_hub_tbl) %>% 
    unnest(cols=c(bins)) %>% 
    mutate(end=bins + res_num[res] - 1)
  
  return(   GenomicRanges::reduce(GRanges(seqnames=tmp_hub_tbl$chr,
                        ranges = IRanges(start=tmp_hub_tbl$bins,
                                         end=tmp_hub_tbl$end
                        ))))
  
})


foot_vec<-unlist(lapply(fdr_cand_Grange_l,function(x){
  sum(width(x))
}))

CAGE_vec<-unlist(lapply(fdr_cand_Grange_l,function(x){
  length(unique(queryHits(findOverlaps(CAGE_GRange,x))))
}))

tibble(FDR=unique(dagger_mres_tbl$FDR),foot=foot_vec,CAGE.n=CAGE_vec) %>% 
  ggplot(.,aes(FDR,foot))+
  geom_line()+
  geom_vline(xintercept = 0.01)+
  theme_classic()
ggsave("~/Documents/Scientia_projects/Files/img/FDR_sensitivity_cl_count.png")
