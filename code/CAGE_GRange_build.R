renv::install("tidyverse")
renv::install("bioc::GenomicRanges")

library(tidyverse)
library(GenomicRanges)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
CAGE_file<-"~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_coord_tbl.Rda"
out_file<-"./data/CAGE_Grange.Rda"
cage_coord_tbl<-get(load(CAGE_file))
cage_coord_tbl<-cage_coord_tbl %>% filter(!(is.na(start)))
cage_chr_Grange<-GRanges(seqnames=cage_coord_tbl$chr,
                            ranges = IRanges(start=as.numeric(cage_coord_tbl$start),
                                             end=as.numeric(cage_coord_tbl$end)
                            ))
save(cage_chr_Grange,file=out_file)
