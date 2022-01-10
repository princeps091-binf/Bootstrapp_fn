#Cage annotation V2
library(vroom)
library(dplyr)
library(tidyr)
library(rtracklayer)
library(readr)
cage<-vroom("~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt",comment = '#')
#Better to only consider the immediate annotation rather than the DPI cluster annotation
cage_ann<-vroom("Documents/multires_bhicect/data/epi_data/DPIcluster_hg19_20120116.permissive_set.GencodeV10_annotated.osc",comment = '#')
tss_ann_b<-tss_ann_b%>%mutate(gene_b=strsplit(tss_ann_b$annotation.GencodeV10.names,split=','))%>%dplyr::select(namess,gene_b)
tss_ann<-cage%>%dplyr::select(`00Annotation`)
colnames(tss_ann)[1]<-"namess"
tss_ann<-tss_ann%>%left_join(.,cage_ann)%>%dplyr::select(namess,annotation.GencodeV10.names)
rm(cage,cage_ann)
save(tss_ann,file = "~/Documents/multires_bhicect/data/epi_data/CAGE/CAGE_tss_ann_smpl.Rda")
load("~/Documents/multires_bhicect/data/epi_data/CAGE/CAGE_tss_ann_smpl.Rda")
tss_ann$tr<-lapply(tss_ann$annotation.GencodeV10.names,function(x)unlist(lapply(strsplit(unlist(strsplit(x,split = ",")),split = "\\."),'[',1)))
unique(grep("^EN",unlist(tss_ann$tr),invert = T,value = T))
# Convert EST to ESG ID
hg19_ENST_to_ENSG <- read_delim("Documents/multires_bhicect/data/epi_data/hg19_ENST_to_ENSG.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)

hg19_ENST_to_ENSG%>%dplyr::select(name,name2)
enst_to_ensg_conv<-hg19_ENST_to_ENSG$name2
names(enst_to_ensg_conv)<-hg19_ENST_to_ENSG$name

cl<-makeCluster(5)
clusterEvalQ(cl, {
  library(dplyr)
})
clusterExport(cl,c('enst_to_ensg_conv'))
tmp_ensg_gene<-parLapply(cl,tss_ann$tr,function(x){
  unique(enst_to_ensg_conv[x])})
stopCluster(cl)
rm(cl)
tmp_ensg_gene<-lapply(tmp_ensg_gene,function(x)x[!(is.na(x))])
tss_ann<-tss_ann%>%mutate(ENSG=tmp_ensg_gene)
tss_ann$ENSG_vec<-unlist(lapply(tss_ann$ENSG,function(x)paste(x,collapse = ',')))
save(tss_ann,file = "~/Documents/multires_bhicect/data/epi_data/CAGE/CAGE_tss_ann_smpl.Rda")
#Compare the annotation composition of CAGE-tss vs full CAGE set for annotations
tss_compo<-tss_ann%>%left_join(.,cage_ann%>%dplyr::select(namess,annotation.GencodeV10.class))%>%group_by(annotation.GencodeV10.class)%>%summarise(n=n())%>%mutate(type="TSS")
all_compo<-cage_ann%>%group_by(annotation.GencodeV10.class)%>%summarise(n=n())%>%mutate(type="ALL")
compo_tbl<-tss_compo%>%bind_rows(.,all_compo)
gg_compo<-ggplot(compo_tbl,aes(type,n,fill=annotation.GencodeV10.class))+geom_bar(stat = "identity",position = "fill")+scale_fill_manual(values = colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(compo_tbl$annotation.GencodeV10.class))))
ggsave("~/Documents/multires_bhicect/weeklies/weekly21/img/cage_tss_ann_compo.png",gg_compo)
