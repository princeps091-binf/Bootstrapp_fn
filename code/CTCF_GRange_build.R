# Convert CTCF bed files to GRange object
library(rtracklayer)
library(GenomicRanges)

CTCF_Grange<-import.bedGraph("./data/ENCFF288RFS.bed.gz")
save(CTCF_Grange,file = "./data/CTCF_Grange.Rda")
