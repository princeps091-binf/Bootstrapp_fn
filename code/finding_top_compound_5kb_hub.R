library(GenomicRanges)
library(tidyverse)
library(furrr)
library(data.tree)
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
#-----------------------------------------
tmp_files<-grep("_5kb_tss_compound_hub.Rda$",list.files("./data/candidate_compound_hub/"),value=T)
tmp_l<-vector("list",length(tmp_files))
names(tmp_l)<-unlist(lapply(strsplit(tmp_files,split="_"),'[',1))
for(file in tmp_files){
  message(strsplit(file,split="_")[[1]][1])
  compound_hub_5kb_tbl<-get_obj_in_fn(paste0("./data/candidate_compound_hub/",file))
  
  top_compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(compound_hub_5kb_tbl$chr),function(chromo){
    tmp_tbl<-compound_hub_5kb_tbl %>% 
      filter(chr==chromo)
    return(tmp_tbl %>% 
             filter(parent.hub %in% tmp_tbl$parent.hub[which(!(tmp_tbl$parent.hub %in% unique(unlist(tmp_tbl$ch.hub))))])
    )
    
  }))
  
  tmp_l[[strsplit(file,split="_")[[1]][1]]]<-  top_compound_hub_5kb_tbl %>% 
    group_by(res) %>% 
    summarise(n=n()) %>% 
    mutate(line=strsplit(file,split="_")[[1]][1])
  
  
}

tmp_tbl<-do.call(bind_rows,tmp_l)
rm(tmp_l)
tmp_tbl %>% 
  mutate(res=fct_relevel(res,names(res_num)),line=fct_relevel(line,c("GM12878","HMEC","H1"))) %>% 
  ggplot(.,aes(line,n,fill=res))+geom_bar(stat="identity",position="fill")+scale_fill_brewer(palette="Set1")
