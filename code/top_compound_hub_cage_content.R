library(GenomicRanges)
library(tidyverse)
library(furrr)
library(data.tree)
library(igraph)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
input_data_fn<-function(tmp_file){
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl) 
}

#-----------------------------------------
compound_hub_5kb_file<-"./data/candidate_compound_hub/GM12878_5kb_tss_compound_hub.Rda"
cage_peak_GRange_file<-"./data/GRanges/CAGE_union_GM12878_Grange.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"

compound_hub_5kb_tbl<-input_data_fn(compound_hub_5kb_file)
cage_GRange<-input_data_fn(cage_peak_GRange_file)
mcols(cage_GRange)<-tibble(ID=paste0("CAGE_",1:length(cage_GRange)))
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
    mutate(bins=map(bins,as.numeric)) %>% 
    mutate(ends=pmap(list(res,bins),function(res,bins){
      return(bins + res_num[res]-1)
    })) %>% 
    mutate(GRange=pmap(list(chr,bins,ends),function(chr,bins,ends){
      return(GRanges(seqnames=chr,
                     ranges = IRanges(start=bins,
                                      end=ends
                     )))
    })) %>% 
    dplyr::select(-c(bins,ends))

}))
plan(multisession, workers=3)

top_compound_hub_5kb_tbl<-top_compound_hub_5kb_tbl %>% 
  mutate(cage.content=future_map(GRange,function(x){
    mcols(cage_GRange)$ID[unique(queryHits(findOverlaps(cage_GRange,x)))]
  }))

top_compound_hub_5kb_tbl %>% 
  group_by(res) %>% 
  summarise(npeak=length(unique(unlist(cage.content))))
#---------------------------------
# Compare with original 5kb compounf set
compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(compound_hub_5kb_tbl$chr),function(chromo){
  message(chromo)
  base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  tmp_tbl<-compound_hub_5kb_tbl %>% 
    filter(chr==chromo) %>% 
    mutate(bins=chr_spec_res$cl_member[parent.hub]) %>% 
    mutate(bins=map(bins,as.numeric)) %>% 
    mutate(ends=pmap(list(res,bins),function(res,bins){
      return(bins + res_num[res]-1)
    })) %>% 
    mutate(GRange=pmap(list(chr,bins,ends),function(chr,bins,ends){
      return(GRanges(seqnames=chr,
                     ranges = IRanges(start=bins,
                                      end=ends
                     )))
    })) %>% 
    dplyr::select(-c(bins,ends))
  
}))
plan(multisession, workers=3)

compound_hub_5kb_tbl<-compound_hub_5kb_tbl %>% 
  mutate(cage.content=future_map(GRange,function(x){
    mcols(cage_GRange)$ID[unique(queryHits(findOverlaps(cage_GRange,x)))]
  }))
