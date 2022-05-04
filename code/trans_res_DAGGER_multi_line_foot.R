library(tidyverse)
library(furrr)
library(data.tree)
library(GenomicRanges)
library(scales)
library(valr)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
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

#-----------------------------------------
# explicit evaluation of ancestry
top_hub_files<-list(H1="./data/DAGGER_tbl/trans_res/H1_union_top_trans_res_dagger_tbl.Rda",
                    GM12878="./data/DAGGER_tbl/trans_res/GM12878_union_top_trans_res_dagger_tbl.Rda",
                    HMEC="./data/DAGGER_tbl/trans_res/HMEC_union_top_trans_res_dagger_tbl.Rda")

spec_res_files<-list(H1="~/Documents/multires_bhicect/data/H1/Dekker/spec_res/",
                     GM12878="~/Documents/multires_bhicect/data/GM12878/spec_res/",
                     HMEC="~/Documents/multires_bhicect/data/HMEC/spec_res/")

#-----------------------------------------
hg19_coord <- read_delim("~/Documents/multires_bhicect/data/hg19.genome", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
names(hg19_coord)<-c("chrom","size")

hub_GRanges_l<-lapply(names(top_hub_files),function(i){
  
  top_hub_tbl<-tbl_in_fn(top_hub_files[[i]]) %>% 
    mutate(res=str_split_fixed(node,"_",2)[,1])
  tmp_res_file<-spec_res_files[i]
  top_hub_tbl<-Build_coord_fn(top_hub_tbl,tmp_res_file) %>% 
    mutate(GRange=pmap(list(chr,res,bins),function(chromo,res,bins){
      Build_GRange_fn(chromo,res,bins,res_num)
    }))
  return(IRanges::reduce(do.call("c",top_hub_tbl$GRange)))
})
names(hub_GRanges_l)<-names(top_hub_files)

tmp_cl_tbl<-hub_GRanges_l[[2]] %>% as_tibble %>% dplyr::select(seqnames,start,end)%>%dplyr::rename(chrom=seqnames)

plan(multisession,workers=5)
rn_hub_foot_tbl<-future_map(1:3,function(x){
  rn_pol<-bed_shuffle(tmp_cl_tbl,genome = hg19_coord,max_tries=1e8,within=T)
  rn_GRange<-GRanges(seqnames=rn_pol$chrom,
                     ranges = IRanges(start=rn_pol$start,
                                      end=rn_pol$end))
  inter_tbl<-tibble(as.data.frame(rn_GRange)) 
  return(inter_tbl %>% mutate(boot=paste0("boot.",x)) %>% dplyr::select(seqnames,start,end,boot))
  
})
plan(sequential)
gg_foot<-do.call(bind_rows,rn_hub_foot_tbl) %>% 
  bind_rows(.,tmp_cl_tbl %>% mutate(boot="obs") %>% dplyr::rename(seqnames=chrom)) %>% 
  mutate(seqnames=fct_relevel(seqnames,paste0('chr',1:22)))%>%
  ggplot(.,aes(xmin=start,xmax=end,ymin=0,ymax=1))+
  geom_rect()+
  facet_grid(seqnames~boot)+
  scale_fill_brewer(palette = "Dark2")

gg_foot<-gg_foot +theme(axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        legend.position="bottom")+
  scale_x_continuous(labels= label_number(scale = 1/1e6,suffix="Mb"))

gg_foot
ggsave("~/Documents/multires_bhicect/weeklies/weekly57/img/top_hub_rn_foot.png",width = 40,height = 23,units = "cm",gg_foot)

hub_coord_tbl<-do.call(bind_rows,lapply(seq_along(hub_GRanges_l),function(x){
  inter_tbl<-tibble(as.data.frame(hub_GRanges_l[[x]])) %>% mutate(line=names(hub_GRanges_l)[x]) 
  return(inter_tbl)
})) 
  
gg_foot<-hub_coord_tbl %>% 
  mutate(seqnames=fct_relevel(seqnames,paste0('chr',1:22)))%>%
  ggplot(.,aes(xmin=start,xmax=end,ymin=0,ymax=1,fill=line))+
  geom_rect()+
  facet_grid(seqnames~line)+
  scale_fill_brewer(palette = "Dark2")

gg_foot<-gg_foot +theme(axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        legend.position="bottom")+
  scale_x_continuous(labels= label_number(scale = 1/1e6,suffix="Mb"))

gg_foot
ggsave("~/Documents/multires_bhicect/weeklies/weekly57/img/top_hub_foot.png",width = 40,height = 23,units = "cm",gg_foot)
