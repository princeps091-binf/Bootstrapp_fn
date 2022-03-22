library(tidyverse)
library(furrr)

#---------------------------------------------------------------------
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
get_obj_in_fn<-function(file){
  out_tbl<-get(load(file))
  tmp_obj<-names(mget(load(file)))
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
      mutate(bins=chr_spec_res$cl_member[parent.hub]) %>% 
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
hub_5kb_file<-"./data/candidate_compound_hub/H1_5kb_tss_compound_hub.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/H1/Dekker/spec_res/"
TAD_file<-"./data/pval_tbl/CAGE_union_H1_TAD_pval_tbl.Rda"

cage_peak_Grange_file<-"./data/GRanges/CAGE_union_H1_Grange.Rda"

TAD_tbl<-get_obj_in_fn(TAD_file)%>% 
  group_by(chr) %>% 
  mutate(FDR=p.adjust(emp.pval,method="fdr")) %>% 
  filter(FDR<=0.01)

cage_GRange<-get_obj_in_fn(cage_peak_Grange_file)
mcols(cage_GRange)<-tibble(ID=paste0("CAGE_",1:length(cage_GRange)))

compound_hub_5kb_tbl<-get_obj_in_fn(hub_5kb_file)

top_compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(compound_hub_5kb_tbl$chr),function(chromo){
  tmp_tbl<-compound_hub_5kb_tbl %>% 
    filter(chr==chromo)
  return(tmp_tbl %>% 
           filter(parent.hub %in% tmp_tbl$parent.hub[which(!(tmp_tbl$parent.hub %in% unique(unlist(tmp_tbl$ch.hub))))])
  )
  
}))
top_compound_hub_5kb_tbl<-Build_coord_fn(top_compound_hub_5kb_tbl,spec_res_file) %>% 
  mutate(GRange=pmap(list(chr,res,bins),function(chromo,res,bins){
    Build_GRange_fn(chromo,res,bins,res_num)
  }))
plan(multisession,workers=3)

top_compound_hub_5kb_tbl<-top_compound_hub_5kb_tbl %>% 
  mutate(peak.content=future_map(GRange,function(x){
    mcols(cage_GRange)$ID[subjectHits(findOverlaps(x,cage_GRange))]
    
  }))
plan(sequential)

plan(multisession,workers=3)

TAD_tbl<-TAD_tbl%>% 
  mutate(peak.content=future_map(GRange,function(x){
    mcols(cage_GRange)$ID[subjectHits(findOverlaps(x,cage_GRange))]
    
  }))
plan(sequential)


TAD_GRange<-reduce(do.call("c",unlist(TAD_tbl$GRange)))
top_hub_GRange<-reduce(do.call("c",unlist(top_compound_hub_5kb_tbl$GRange)))

sum(width(intersect(TAD_GRange,top_hub_GRange)))/sum(width(reduce(c(TAD_GRange,top_hub_GRange))))
sum(width(intersect(TAD_GRange,top_hub_GRange)))/sum(width(top_hub_GRange))
sum(width(intersect(TAD_GRange,top_hub_GRange)))/sum(width(TAD_GRange))


up_list<-list(TAD=unique(unlist(TAD_tbl$peak.content)),top.hub=unique(unlist(top_compound_hub_5kb_tbl$peak.content)))

library(UpSetR)
upset(fromList(up_list),order.by = "freq")


