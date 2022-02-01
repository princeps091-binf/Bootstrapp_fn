# Empirical p-value process re-factored as functions
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
feature_annotation_fn<-function(txdb,chr_feature_Grange,fn_file){
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
  rm(ChIPseekerEnv,envir = globalenv())
  rm(n5,n3,nexon,nintron,n1kb,n2kb,n3kb,ndown,ninter)
  n_vec<-n_vec[n_vec>0]
  return(n_vec)
  
}

rn_feature_GRange_build_fn<-function(n_vec,fn_bed_l,hg19_coord,tmp_cage_tbl,fn_file){
  
  fn_env<-environment()
  
  rn_fn_coord_l<-vector('list',length(n_vec))
  names(rn_fn_coord_l)<-names(n_vec)
  for(f in names(n_vec)){
    # 
    tmp_n<-n_vec[f]
    if(f==fn_file[5]){
      cl<-makeCluster(5)
      clusterEvalQ(cl, {
        library(dplyr)
        library(valr)
        print("node ready")
      })
      clusterExport(cl,c("f","tmp_n","fn_bed_l","tmp_cage_tbl","hg19_coord"),envir = fn_env)
      
      rn_fn_coord_l[[f]]<-parLapply(cl,1:100,function(x){
        rn_pol<-bed_shuffle(tmp_cage_tbl%>%sample_n(tmp_n),genome = hg19_coord,excl = fn_bed_l[[f]],within = T,max_tries=1e6)
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
      clusterExport(cl,c("f","tmp_n","fn_bed_l","tmp_cage_tbl","hg19_coord"),envir = fn_env)
      rn_fn_coord_l[[f]]<-parLapply(cl,1:100,function(x){
        # Try pattern to not abort at the shuffling step
        ## Sampleing prior to shuffling ensures a greater likelihood of success
        rn_pol<-try(valr::bed_shuffle(x = tmp_cage_tbl%>%sample_n(tmp_n),genome = hg19_coord,incl = fn_bed_l[[f]],within=T,max_tries=1e6),silent=T)
        return(rn_pol)
      })
      stopCluster(cl)
      rm(cl)
      # Collect the successful shuffling by eliminating the shuffles producing try-error objects
      rn_fn_coord_l[[f]]<-rn_fn_coord_l[[f]][!(unlist(lapply(rn_fn_coord_l[[f]],function(x)any(class(x) %in% "try-error"))))]
    }
  }
  
  #Assemble these blocks into 10000 rnadom combo
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(dplyr)
    print("node ready")
  })
  clusterExport(cl,c("rn_fn_coord_l"),envir = fn_env)
  
  rn_peak_coord_tbl_l<-parLapply(cl,1:10000,function(x){
    return(do.call(bind_rows,lapply(rn_fn_coord_l,function(f)f[[sample(1:length(f),1)]])))
  })
  stopCluster(cl)
  rm(cl)
  rm(rn_fn_coord_l,envir = fn_env)
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
  rm(rn_peak_coord_tbl_l,envir = fn_env)
  
  return(rn_Grange_l)
  
}

empirical_pval_compute_fn<-function(chromo,cl_folder,cl_file,feature_Grange,fn_repo,txdb,hg19_coord){

  main_fn_env<-environment()
  
  
  chr_feature_Grange<-feature_Grange[seqnames(feature_Grange)==chromo]
  
  cat(green(chromo), " Bootstrapping started \n")
  
  fn_folder<-paste0(fn_repo,chromo,"/")
  fn_file<-grep('BED$',list.files(fn_folder),value = T)
  fn_bed_l<-lapply(fn_file,function(f){
    read_bed(paste0(fn_folder,f),n_fields = 3)
  })
  names(fn_bed_l)<-fn_file
  txdb_chr <- txdb
  seqlevels(txdb_chr)  <- chromo
  
  cat(green(chromo)," feature annotation \n")
  
  n_vec<-feature_annotation_fn(txdb_chr,chr_feature_Grange,fn_file)
  
  cl_chr_tbl<-get(load(paste0(cl_folder,chromo,cl_file)))
  tmp_obj<-names(mget(load(paste0(cl_folder,chromo,cl_file))))
  rm(list=tmp_obj)
  rm(tmp_obj)
  
  cat(green(chromo)," build GRangeList object \n")
  
  #table collecting the observed CAGE-peak coordinates
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(GenomicRanges)
    print("node ready")
  })
  clusterExport(cl,c("chr_feature_Grange"),envir=main_fn_env)#
  cl_inter_vec<-unlist(parLapply(cl,cl_chr_tbl$GRange,function(x){
    sum(countOverlaps(x,chr_feature_Grange))
  }))
  stopCluster(cl)
  rm(cl)
  
  cl_chr_tbl<-cl_chr_tbl%>%mutate(feature_n=cl_inter_vec)%>% filter(feature_n>0)
  
  cl_list<-GRangesList(cl_chr_tbl$GRange)  
  
  tmp_cage_tbl<-chr_feature_Grange %>% as_tibble %>% dplyr::select(seqnames,start,end)%>%dplyr::rename(chrom=seqnames)
  
  cat(green(chromo)," build random feature coordinates \n")
  
  rn_Grange_l<-rn_feature_GRange_build_fn(n_vec,fn_bed_l,hg19_coord,tmp_cage_tbl,fn_file)
  
  cat(green(chromo)," compute empirical p-value \n")
  
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(GenomicRanges)
    print("node ready")
  })
  clusterExport(cl,c("cl_list"),envir=main_fn_env)
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
  
  return(cl_chr_tbl%>%mutate(emp.pval=cl_emp_pval))
  
}

#-----------------------------------------
cl_folder<-"./data/GRanges/BHiCect_Grange/GM12878/"
cl_file<-"_BHiCect_cl.Rda"
feature_file<-"./data/GRanges/CAGE_union_GM12878_Grange.Rda"
out_file<-"./data/pval_tbl/CAGE_union_GM12878_pval_tbl.Rda"

feature_Grange<-get(load(feature_file))
tmp_obj<-names(mget(load(feature_file)))
rm(list=tmp_obj)
rm(tmp_obj)

hg19_coord <- read_delim("~/Documents/multires_bhicect/data/hg19.genome", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
names(hg19_coord)<-c("chrom","size")

fn_repo<-"~/Documents/multires_bhicect/data/epi_data/fn_BED/"

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

chr_set<-list.files(fn_repo)
cl_chr_emp_pval_l<-lapply(chr_set,function(chromo){
  cl_chr_tbl<-empirical_pval_compute_fn(chromo,cl_folder,cl_file,feature_Grange,fn_repo,txdb,hg19_coord)
  seqlevels(txdb)<-seqlevels0(txdb)
  return(cl_chr_tbl)
})

cl_chr_emp_pval_tbl<-do.call(bind_rows,cl_chr_emp_pval_l)

save(cl_chr_emp_pval_tbl,file=out_file)
