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
#-------------------------------------------------------------------------------------------------------
compound_hub_5kb_file<-"./data/candidate_compound_hub/H1_5kb_tss_compound_hub.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/H1/Dekker/spec_res/"
gene_GRange_file<-"~/Documents/multires_bhicect/GO_enrichment_viz/data/CAGE_H1_gene_GRange.Rda"
gene_conv_tbl_file<-"~/Documents/multires_bhicect/GO_enrichment_viz/data/gene_name_conv_tbl.Rda"
semantic_module_file<-"~/Documents/multires_bhicect/GO_enrichment_viz/data/semantic_module_tbl/H1_semantic_module_tbl.Rda"

gene_GRange<-tbl_in_fn(gene_GRange_file)

compound_hub_5kb_tbl<-tbl_in_fn(compound_hub_5kb_file)

gene_conv_tbl<-tbl_in_fn(gene_conv_tbl_file)

sem_mod_tbl<-tbl_in_fn(semantic_module_file)
#-------------------------------------------------------------------------------------------------------
## Subset the mutually exclusive collection of compound hubs
top_compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(compound_hub_5kb_tbl$chr),function(chromo){
  tmp_tbl<-compound_hub_5kb_tbl %>% 
    filter(chr==chromo)
  return(tmp_tbl %>% 
           filter(parent.hub %in% tmp_tbl$parent.hub[which(!(tmp_tbl$parent.hub %in% unique(unlist(tmp_tbl$ch.hub))))])
  )
  
}))
## Build GRange object for each top hub
top_compound_hub_5kb_tbl<-Build_coord_fn(top_compound_hub_5kb_tbl,spec_res_file) %>% 
    mutate(GRange=pmap(list(chr,res,bins),function(chromo,res,bins){
      Build_GRange_fn(chromo,res,bins,res_num)
    }))
## Compute the ENSG content for each top hub
top_compound_hub_5kb_tbl<-top_compound_hub_5kb_tbl %>% 
  mutate(ENSG.content=map(GRange,function(x){
    unique(unlist(mcols(gene_GRange)$ENSG[queryHits(findOverlaps(gene_GRange,x))]))
  }))
## Convert the ENSG content to entrez gene content
top_compound_hub_5kb_tbl<-top_compound_hub_5kb_tbl %>% 
  mutate(entrez.content=map(ENSG.content,function(x){
    return(gene_conv_tbl %>% 
             filter(ensembl_gene_id %in% x) %>% 
             dplyr::select(entrezgene_id) %>% 
             distinct %>% unlist)
  }))

mod_n<-3
sem_mod_tbl$GO.tbl[[mod_n]] %>% arrange(FDR)

tmp_mod_entrez_set<-sem_mod_tbl$entrez.content[[mod_n]]


top_compound_hub_5kb_tbl<-top_compound_hub_5kb_tbl %>% 
  mutate(mod.content=map_int(entrez.content,function(x){
    sum(tmp_mod_entrez_set %in% as.character(x))
  })) 

top_compound_hub_5kb_tbl %>% dplyr::select(chr,parent.hub,mod.content) %>% arrange(desc(mod.content))

top_compound_hub_5kb_coord_tbl<-top_compound_hub_5kb_tbl %>% 
  mutate(coord=pmap(list(parent.hub,GRange,mod.content),function(parent.hub,GRange,mod.content){
    tibble(as.data.frame(GRange)) %>% 
#      mutate(hub=parent.hub,content=mod.content/length(tmp_mod_entrez_set)) %>% 
      mutate(hub=parent.hub,content=ifelse(mod.content>0,"in","out"))
  })) %>% 
  dplyr::select(chr,parent.hub,coord)
coord_tbl<-do.call(bind_rows,top_compound_hub_5kb_coord_tbl$coord)

gg_foot<-coord_tbl%>%
  filter(grepl("^5kb_",hub)) %>% 
  mutate(seqnames=fct_relevel(seqnames,paste0('chr',1:22)))%>%
  ggplot(.,aes(xmin=start,xmax=end,ymin=0,ymax=1,fill=content))+
  geom_rect()+
  facet_grid(seqnames~.)+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()


gg_foot<-gg_foot +theme(axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank())+
  scale_x_continuous(labels= label_number(scale = 1/1e6,suffix="Mb"))

gg_foot

coord_tbl %>% 
  filter(grepl("^5kb_",hub)) %>% 
  group_by(seqnames,hub) %>% 
  summarise(w=sum(width)) %>% 
  ggplot(.,aes(w))+geom_density()+scale_x_log10()
#----------------
mod_n<-3
sem_mod_tbl$GO.tbl[[mod_n]] %>% arrange(FDR)

tmp_mod_entrez_set<-sem_mod_tbl$entrez.content[[mod_n]]

tmp_chi_tbl<-top_compound_hub_5kb_tbl %>% 
  mutate(mod.content=map_int(entrez.content,function(x){
    sum(tmp_mod_entrez_set %in% as.character(x))
  })) %>% filter(mod.content>0) %>% 
  group_by(res) %>% 
  summarise(mod.n=sum(mod.content),entrez.n=length(unique(unlist(entrez.content))))
tmp_obs<-tmp_chi_tbl$mod.n
names(tmp_obs)<-tmp_chi_tbl$res
exp_p<-tmp_chi_tbl$entrez.n/sum(tmp_chi_tbl$entrez.n)
chi_obj<-chisq.test(tmp_obs,p=exp_p)
chi_obj$p.value
chi_obj$observed - chi_obj$expected
#-----------------------------------

up_l<-lapply(sem_mod_tbl$sem.mod,function(mod_n){

  tmp_mod_entrez_set<-sem_mod_tbl$entrez.content[[mod_n]]
  
  
  top_compound_hub_5kb_tbl %>% 
    mutate(mod.content=map_int(entrez.content,function(x){
      sum(tmp_mod_entrez_set %in% as.character(x))
    })) %>%
    filter(mod.content>0) %>% 
    mutate(sem_mod=mod_n,ID=paste(chr,parent.hub,sep="_")) %>% 
    dplyr::select(ID) %>% unlist
  
})
names(up_l)<-paste0("sem.mod.",1:length(up_l))
library(UpSetR)
upset(fromList(up_l),order.by = "freq",nsets = length(up_l))
inter_tbl<-as_tibble(fromList(up_l))
inter_tbl<-inter_tbl %>% 
  mutate(hubs=unique(unlist(up_l)))
inter_tbl %>% 
#  filter(sem.mod.1 == 1 & sem.mod.2 == 1 &sem.mod.3 == 1 &sem.mod.4 == 1 &sem.mod.5 == 1 &sem.mod.6 == 1 ) %>% 
  filter(sem.mod.1 == 1 & sem.mod.2 == 1 &sem.mod.3 == 1 &sem.mod.4 == 1  ) %>% 
  mutate(res=map_chr(hubs,function(x)strsplit(x,split="_")[[1]][2])) %>% 
  group_by(res) %>% 
  summarise(n=n())

top_compound_hub_5kb_tbl %>% 
  group_by(res) %>% 
  summarise(n=n())

inter_tbl %>% 
  mutate(chr=str_split_fixed(hubs,pattern = "_",n=2)[,1]) %>% 
  mutate(parent.hub=str_split_fixed(hubs,pattern = "_",n=2)[,2]) %>% 
  dplyr::select(-hubs) %>% 
  #  filter(sem.mod.1 == 1 & sem.mod.2 == 1 &sem.mod.3 == 1 &sem.mod.4 == 1 &sem.mod.5 == 1 &sem.mod.6 == 1 ) %>% 
  filter(sem.mod.1 == 1 & sem.mod.2 == 1 &sem.mod.3 == 1 &sem.mod.4 == 1  ) %>% 
  dplyr::select(chr,parent.hub) %>% 
  mutate(set="multi-function") %>% 
  full_join(.,top_compound_hub_5kb_tbl %>% 
              dplyr::select(chr,parent.hub,peak.content)) %>% 
  mutate(set=ifelse(is.na(set),'out',set)) %>% 
  mutate(res=str_split_fixed(parent.hub,pattern = "_",n=4)[,1]) %>% 
  dplyr::select(parent.hub,res,set) %>% 
  distinct() %>% 
  group_by(res,set) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(x=res,n,fill=set))+
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()


inter_tbl %>% 
  mutate(chr=str_split_fixed(hubs,pattern = "_",n=2)[,1]) %>% 
  mutate(parent.hub=str_split_fixed(hubs,pattern = "_",n=2)[,2]) %>% 
  dplyr::select(-hubs) %>% 
#  filter(sem.mod.1 == 1 & sem.mod.2 == 1 &sem.mod.3 == 1 &sem.mod.4 == 1 &sem.mod.5 == 1 &sem.mod.6 == 1 ) %>% 
  filter(sem.mod.1 == 1 & sem.mod.2 == 1 &sem.mod.3 == 1 &sem.mod.4 == 1  ) %>% 
  dplyr::select(chr,parent.hub) %>% 
  mutate(set="multi-function") %>% 
  full_join(.,top_compound_hub_5kb_tbl %>% 
              dplyr::select(chr,parent.hub,peak.content)) %>% 
  mutate(set=ifelse(is.na(set),'out',set)) %>% 
  mutate(res=str_split_fixed(parent.hub,pattern = "_",n=4)[,1]) %>% 
  dplyr::select(peak.content,res,set) %>% 
  unnest(cols=c(peak.content)) %>% 
  distinct() %>% 
  group_by(res,set) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(x=res,n,fill=set))+
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()
  
