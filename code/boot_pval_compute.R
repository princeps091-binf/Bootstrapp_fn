
library(tidyverse)
library(GenomicRanges)
library(valr)
library(parallel)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(crayon)

options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
## Utility functions

#-----------------------------------------
## Files with cluster, feature GRange objects
cl_folder<-"./data/GRanges/BHiCect_Grange/HMEC/"
cl_file<-"_BHiCect_cl.Rda"
feature_file<-"./data/GRanges/CAGE_enh_HMEC_Grange.Rda"
out_file<-"~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_enh_emp_pval_tbl.Rda"
## GRange for feature of interest
feature_Grange<-get(load(feature_file))
tmp_obj<-names(mget(load(feature_file)))
rm(list=tmp_obj)
rm(tmp_obj)

hg19_coord <- read_delim("~/Documents/multires_bhicect/data/hg19.genome", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
names(hg19_coord)<-c("chrom","size")

fn_repo<-"~/Documents/multires_bhicect/data/epi_data/fn_BED/"
chr_set<-list.files(fn_repo)

chr_res_l<-vector("list",length(chr_set))
names(chr_res_l)<-chr_set

for (chromo in chr_set){
  message(green(chromo))
  chr_feature_Grange<-feature_Grange[seqnames(feature_Grange)==chromo]
  #Generate random CAGE coordinate
  fn_folder<-paste0(fn_repo,chromo,"/")
  fn_file<-grep('BED$',list.files(fn_folder),value = T)
  fn_bed_l<-lapply(fn_file,function(f){
    read_bed(paste0(fn_folder,f),n_fields = 3)
  })
  names(fn_bed_l)<-fn_file
  
  #Generate the random peak coordinates
  message(green(chromo)," feature annotation")
  
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(ChIPseeker)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  seqlevels(txdb)  <- chromo
  peakAnno <- annotatePeak(chr_feature_Grange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)
  
  rn_annotation<-sample(peakAnno@annoStat$Feature,size = length(chr_feature_Grange),prob = peakAnno@annoStat$Frequency/100,replace = T)
  #check number of peaks from that category
  n5<-length(grep("5'",as.character(rn_annotation)))
  n3<-length(grep("3'",as.character(rn_annotation)))
  nexon<-length(grep("Exon",as.character(rn_annotation)))
  nintron<-length(grep("Intron",as.character(rn_annotation)))
  n1kb<-length(grep("1kb",as.character(rn_annotation)))
  n2kb<-length(grep("2kb",as.character(rn_annotation)))
  n3kb<-length(grep("3kb",as.character(rn_annotation)))
  ndown<-length(grep("Down",as.character(rn_annotation)))
  ninter<-length(grep("Inter",as.character(rn_annotation)))
  n_vec<-c(n3,n5,ndown,nexon,ninter,nintron,n1kb,n2kb,n3kb)
  names(n_vec)<-fn_file
  rm(ChIPseekerEnv)
  rm(n5,n3,nexon,nintron,n1kb,n2kb,n3kb,ndown,ninter)
  n_vec<-n_vec[n_vec>0]
  #----------------------------------------------------
  cl_chr_tbl<-get(load(paste0(cl_folder,chromo,cl_file)))
  tmp_obj<-names(mget(load(paste0(cl_folder,chromo,cl_file))))
  rm(list=tmp_obj)
  rm(tmp_obj)
  
  
  #cl_chr_tbl<-cl_Grange_tbl %>% filter(chr==chromo)
  message(green(chromo)," build GRangeList object")
  
  cl_list<-GRangesList(cl_chr_tbl$GRange)  
  #table collecting the observed CAGE-peak coordinates
  cl_chr_tbl<-cl_chr_tbl%>%mutate(feature_n=countOverlaps(cl_list,chr_feature_Grange))%>% filter(feature_n>0)
  
  tmp_cage_tbl<-chr_feature_Grange %>% as_tibble %>% dplyr::select(seqnames,start,end)%>%dplyr::rename(chrom=seqnames)
  #Build random coord sets for each categories
  message(green(chromo)," build random feature coordinates")
  
  rn_fn_coord_l<-vector('list',length(n_vec))
  names(rn_fn_coord_l)<-names(n_vec)
  for(f in names(n_vec)){
    cat(f)
    tmp_n<-n_vec[f]
    if(f==fn_file[5]){
      cl<-makeCluster(5)
      clusterEvalQ(cl, {
        library(dplyr)
        library(valr)
        print("node ready")
      })
      clusterExport(cl,c("f","tmp_n","fn_bed_l","tmp_cage_tbl","hg19_coord"))
      
      rn_fn_coord_l[[f]]<-parLapply(cl,1:100,function(x){
        rn_pol<-bed_shuffle(tmp_cage_tbl,genome = hg19_coord,excl = fn_bed_l[[f]],within = T,max_tries=1e3)%>%sample_n(tmp_n)
        return(rn_pol)
      })
      stopCluster(cl)
      rm(cl)
      
    }
    if(f!=fn_file[5]){
      cl<-makeCluster(5)
      clusterEvalQ(cl, {
        library(dplyr)
        library(valr)
        print("node ready")
      })
      clusterExport(cl,c("f","tmp_n","fn_bed_l","tmp_cage_tbl","hg19_coord"))
      rn_fn_coord_l[[f]]<-parLapply(cl,1:100,function(x){
        rn_pol<-valr::bed_shuffle(x = tmp_cage_tbl,genome = hg19_coord,incl = fn_bed_l[[f]],within=T,max_tries=1e3)%>%sample_n(tmp_n)
        return(rn_pol)
      })
      stopCluster(cl)
      rm(cl)
    }
  }
  
  #Assemble these blocks into 10000 rnadom combo
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(dplyr)
    print("node ready")
  })
  clusterExport(cl,c("rn_fn_coord_l"))
  
  rn_peak_coord_tbl_l<-parLapply(cl,1:10000,function(x){
    return(do.call(bind_rows,lapply(rn_fn_coord_l,function(f)f[[sample(1:length(f),1)]])))
  })
  stopCluster(cl)
  rm(cl)
  
  #Convert to Grange object
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(GenomicRanges)
    library(dplyr)
    print("node ready")
  })
  rn_Grange_l<-parLapply(cl,rn_peak_coord_tbl_l,function(x){
    rnp_Grange<-GRanges(seqnames=x$chrom,
                        ranges = IRanges(start=x$start,
                                         end=x$end)
    )
    return(rnp_Grange)
  })
  stopCluster(cl)
  rm(cl)
  rm(rn_peak_coord_tbl_l)
  #----------------------------------------------------
  #Compute intersection with observed clusters
  message(green(chromo)," compute empirical p-value")
  
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(GenomicRanges)
    print("node ready")
  })
  clusterExport(cl,c("cl_list"))
  rn_pval_l<-parLapply(cl,rn_Grange_l,function(x){
    countOverlaps(cl_list,x)
  })
  stopCluster(cl)
  rm(cl)
  rm(rn_Grange_l)
  rn_count_mat<-do.call(cbind,rn_pval_l)
  rm(rn_pval_l)
  cl_emp_pval<-unlist(lapply(1:nrow(cl_chr_tbl),function(x){
    (sum(rn_count_mat[x,]>=cl_chr_tbl$feature_n[x])+1)/(ncol(rn_count_mat)+1)
  }))
  
  chr_res_l[[chromo]]<-cl_chr_tbl%>%mutate(emp.pval=cl_emp_pval)
  rm(rn_count_mat,cl_emp_pval) 
  
  
}

cl_emp_pval_tbl<-do.call(bind_rows,chr_res_l)
#cl_emp_pval_tbl<-cl_emp_pval_tbl %>% mutate(FDR=p.adjust(emp.pval,method = 'fdr'))
save(cl_emp_pval_tbl,file=out_file)
