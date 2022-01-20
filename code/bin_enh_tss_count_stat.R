# Examine bin specific distribution of enh and tss
library(tidyverse)
library(GenomicRanges)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
##input hic data of choice
hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat <- read_delim(paste0(dat_file,cl_res,"/",chromo,".txt"), 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.nan(X3)))%>%filter(X1!=X2))
}

Build_GRange_fn<-function(bin_set,res,chr,res_num){
  
  return(GRanges(seqnames=chr,
                 ranges = IRanges(start=as.numeric(bin_set),
                                  end=as.numeric(bin_set)+res_num[res]-1
                 )))
  
}

chr_bin_GRange_fn<-function(dat_file,chromo,cl_res,res_num){
  chr_dat<-hic_dat_in(dat_file,cl_res,chromo)
  bin_set<-unique(c(chr_dat$X1,chr_dat$X2))
  return(Build_GRange_fn(bin_set,cl_res,chromo,res_num))
}
#-----------------------------------------
enh_Grange_file<-"./data/GRanges/CAGE_enh_HMEC_Grange.Rda"
tss_Grange_file<-"./data/GRanges/CAGE_tss_HMEC_Grange.Rda"

enh_Grange<-get(load(enh_Grange_file))
tmp_obj<-names(mget(load(enh_Grange_file)))
rm(list=tmp_obj)
rm(tmp_obj)

tss_Grange<-get(load(tss_Grange_file))
tmp_obj<-names(mget(load(tss_Grange_file)))
rm(list=tmp_obj)
rm(tmp_obj)

dat_file<-"~/Documents/multires_bhicect/data/HMEC/"

cl_res<-"1Mb"
chr_set<-unlist(lapply(strsplit(grep("chr",list.files(paste0(dat_file,"/",cl_res)),value=T),split="\\."),'[',1))

res_tbl<-do.call(bind_rows,lapply(chr_set,function(chromo){
  message(chromo)
  chr_bin_GRange<-chr_bin_GRange_fn(dat_file,chromo,cl_res,res_num)
  
  bin_count_tbl<-tibble(chr=chromo,res=cl_res,bin=start(chr_bin_GRange),enh.n=countOverlaps(chr_bin_GRange,enh_Grange),tss.n=countOverlaps(chr_bin_GRange,tss_Grange))
  
}))

res_tbl %>% ggplot(.,aes(enh.n,tss.n))+geom_point()+scale_y_log10()
