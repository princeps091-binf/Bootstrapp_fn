#Generate the BED files for peak annotation
library(readr)
library(Matrix)
library(data.tree)
library(ggplot2)
library(viridis)
library(rtracklayer)
library(MASS)
library(dplyr)
library(tidyr)
library(caret)
library(bedtoolsr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#-----------------------------------------------------------------------------------------------------------
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------------------------------------------------------------------------
chr_set<-paste0("chr",1:22)
for(i in chr_set){
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  print(paste("producing Functional GRanges for",i))
  seqlevels(txdb) <- i
  chr_prom<-promoters(txdb,upstream = 1000,downstream = 0)
  chr_prom_1_2<-flank(chr_prom,start = T,width = 1000)
  chr_prom_2_3<-flank(chr_prom_1_2,start = T,width = 1000)
  chr_down<-flank(genes(txdb),start=F,width = 1000)
  #exon
  chr_exon<-exons(txdb)
  #intron
  chr_intron<-intronicParts(txdb)
  #5'UTR
  chr_5utr<-unique(unlist(fiveUTRsByTranscript(txdb)))
  #3'UTR
  chr_3utr<-unique(unlist(threeUTRsByTranscript(txdb)))
  #Complement of all these (exclude centroemeres, telomeres and already generated peaks)
  chr_inter_no<-c(chr_5utr,chr_3utr,chr_down,granges(chr_intron),granges(chr_exon),granges(chr_prom_2_3),granges(chr_prom_1_2),granges(chr_prom))
  #exclude the centromeres and telomeres
  gap<-read_delim("~/Documents/multires_bhicect/data/gap.bed",delim = "\t",col_names = F)
  gap<-gap%>%filter(X2==i)%>%arrange(X3)
  gap[1,3]<-1
  gap_Grange<- GRanges(seqnames=gap$X2,
                       ranges = IRanges(start=gap$X3,
                                        end=gap$X4))
  chr_inter_no<-c(chr_inter_no,gap_Grange)
  print(paste("saving BED files for",i))
  #Create BED files for each category
  export(chr_prom,con = paste0("~/Documents/multires_bhicect/data/epi_data/fn_BED/",i,"/",i,"_prom.BED"),format = "bedgraph")
  export(chr_prom_1_2,con = paste0("~/Documents/multires_bhicect/data/epi_data/fn_BED/",i,"/",i,"_prom12.BED"),format = "bedgraph")
  export(chr_prom_2_3,con = paste0("~/Documents/multires_bhicect/data/epi_data/fn_BED/",i,"/",i,"_prom23.BED"),format = "bedgraph")
  export(chr_exon,con = paste0("~/Documents/multires_bhicect/data/epi_data/fn_BED/",i,"/",i,"_exon.BED"),format = "bedgraph")
  export(chr_intron,con = paste0("~/Documents/multires_bhicect/data/epi_data/fn_BED/",i,"/",i,"_intron.BED"),format = "bedgraph")
  export(chr_down,con = paste0("~/Documents/multires_bhicect/data/epi_data/fn_BED/",i,"/",i,"_down.BED"),format = "bedgraph")
  export(chr_5utr,con = paste0("~/Documents/multires_bhicect/data/epi_data/fn_BED/",i,"/",i,"_5utr.BED"),format = "bedgraph")
  export(chr_3utr,con = paste0("~/Documents/multires_bhicect/data/epi_data/fn_BED/",i,"/",i,"_3utr.BED"),format = "bedgraph")
  export(chr_inter_no,con = paste0("~/Documents/multires_bhicect/data/epi_data/fn_BED/",i,"/",i,"_inter_no.BED"),format = "bedgraph")

}

epi_folder<-'~/Documents/multires_bhicect/data/epi_data/H1/'
bed_file<-grep('bed$',list.files(epi_folder),value = T)
load(paste0(epi_folder,'max_POL2_chr_coord_epi.Rda'))
pol2_bed<-import(con = paste0(epi_folder,bed_file[2]),format = "bedgraph")

for(i in chr_set){
  chr_pol2_bed<-pol2_bed[seqnames(pol2_bed)==i]
  
  print(paste("producing Pol2 GRanges for",i))
  export(chr_pol2_bed,con = paste0("~/Documents/multires_bhicect/data/epi_data/H1/POL2_BED/",i,"_pol2_peak.BED"),format = "bedgraph")
  
}
