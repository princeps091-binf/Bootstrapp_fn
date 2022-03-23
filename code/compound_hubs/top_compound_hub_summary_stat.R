# Top-hub generic topo stat
library(tidyverse)
library(GenomicRanges)
library(scales)
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

topo_summary_fn<-function(top_compound_hub_5kb_tbl,spec_res_file,cage_GRange){
  message("Build GRange")
  top_compound_hub_5kb_tbl<-Build_coord_fn(top_compound_hub_5kb_tbl,spec_res_file) %>% 
    mutate(GRange=pmap(list(chr,res,bins),function(chromo,res,bins){
      Build_GRange_fn(chromo,res,bins,res_num)
    }))
  message("Compute Peak content")
  
  ## Compute the ENSG content for each top hub
  top_compound_hub_5kb_tbl<-top_compound_hub_5kb_tbl %>% 
    mutate(peak.content=map_int(GRange,function(x){
      length(unique(queryHits(findOverlaps(cage_GRange,x))))
    }))
  message("Compute topo stats")
  
  tmp_tbl<-top_compound_hub_5kb_tbl %>% 
    mutate(foot=map_int(GRange,function(x)sum(width(x)))) %>%
    mutate(chunk.n=map_int(GRange,function(x)length(x))) %>% 
    mutate(res=fct_relevel(res,names(res_num))) %>% 
    dplyr::select(chr,res,parent.hub,peak.content,foot,chunk.n)
  return(tmp_tbl)
  
}

#-------------------------------------------------------------------------------------------------------
compound_hub_5kb_files<-list(H1="./data/candidate_compound_hub/H1_5kb_tss_compound_hub.Rda",
                             GM12878="./data/candidate_compound_hub/GM12878_5kb_tss_compound_hub.Rda",
                             HMEC="./data/candidate_compound_hub/HMEC_5kb_tss_compound_hub.Rda")

spec_res_files<-list(H1="~/Documents/multires_bhicect/data/H1/Dekker/spec_res/",
                     GM12878="~/Documents/multires_bhicect/data/GM12878/spec_res/",
                     HMEC="~/Documents/multires_bhicect/data/HMEC/spec_res/")

cage_peak_Grange_files<-list(H1="./data/GRanges/CAGE_union_H1_Grange.Rda",
                             GM12878="./data/GRanges/CAGE_union_GM12878_Grange.Rda",
                             HMEC="./data/GRanges/CAGE_union_HMEC_Grange.Rda")

#-------------------------------------------------------------------------------------------------------
## Subset the mutually exclusive collection of compound hubs
topo_summary_tbl<-do.call(bind_rows,lapply(1:length(compound_hub_5kb_files),function(cell_line){
  cage_GRange<-tbl_in_fn(cage_peak_Grange_files[[cell_line]])
  
  compound_hub_5kb_tbl<-tbl_in_fn(compound_hub_5kb_files[[cell_line]])
  
  top_compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(compound_hub_5kb_tbl$chr),function(chromo){
    tmp_tbl<-compound_hub_5kb_tbl %>% 
      filter(chr==chromo)
    return(tmp_tbl %>% 
             filter(parent.hub %in% tmp_tbl$parent.hub[which(!(tmp_tbl$parent.hub %in% unique(unlist(tmp_tbl$ch.hub))))])
    )
    
  }))
  
  tmp_tbl<-topo_summary_fn(top_compound_hub_5kb_tbl,spec_res_files[[cell_line]],cage_GRange) %>% 
    mutate(cell.line=names(spec_res_files)[cell_line])
  return(tmp_tbl)
}))


topo_summary_tbl %>% 
  group_by(cell.line,res) %>% 
  summarise(n=n(),foot=sum(foot),peak=sum(peak.content)) %>% 
  mutate(res=fct_relevel(res,names(res_num))) %>% 
  ggplot(.,aes(cell.line,n,fill=res))+
  geom_bar(stat="identity",position="fill")+
  ylab("number of clusters")+
  scale_fill_viridis_d()
ggsave(filename = "~/Documents/multires_bhicect/weeklies/weekly53/img/top_hub_ncl_bar.png")

topo_summary_tbl %>% 
  group_by(cell.line,res) %>% 
  summarise(n=n(),foot=sum(foot),peak=sum(peak.content)) %>% 
  mutate(res=fct_relevel(res,names(res_num))) %>% 
  ggplot(.,aes(cell.line,foot,fill=res))+
  geom_bar(stat="identity",position="fill")+
  ylab("cluster footprint")+
  scale_fill_viridis_d()
ggsave(filename = "~/Documents/multires_bhicect/weeklies/weekly53/img/top_hub_foot_bar.png")


topo_summary_tbl %>% 
  group_by(cell.line,res) %>% 
  summarise(n=n(),foot=sum(foot),peak=sum(peak.content)) %>% 
  mutate(res=fct_relevel(res,names(res_num))) %>% 
  ggplot(.,aes(cell.line,peak,fill=res))+
  ylab("CAGE-peak content")+
  geom_bar(stat="identity",position="fill")+
  scale_fill_viridis_d()
ggsave(filename = "~/Documents/multires_bhicect/weeklies/weekly53/img/top_hub_npeak_bar.png")


topo_summary_tbl %>% 
  ggplot(.,aes(foot,color=cell.line))+
  geom_density()+
  scale_x_log10()+
  facet_wrap(res~.,scales="free")+
  scale_color_brewer(palette="Set1")
ggsave(filename = "~/Documents/multires_bhicect/weeklies/weekly53/img/top_hub_foot_density.png")

topo_summary_tbl %>% 
  ggplot(.,aes(peak.content,color=cell.line))+
  geom_density()+
  facet_wrap(res~.,scales="free")+
  scale_color_brewer(palette="Set1")
ggsave(filename = "~/Documents/multires_bhicect/weeklies/weekly53/img/top_hub_peak_content_density.png")

topo_summary_tbl %>% 
  ggplot(.,aes(chunk.n,color=cell.line))+
  geom_density()+
  facet_wrap(res~.,scales="free")+
  scale_color_brewer(palette="Set1")

