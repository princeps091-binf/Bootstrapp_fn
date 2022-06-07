library(tidyverse)
library(Matrix)
library(furrr)
library(viridis)
library(data.tree)
library(vroom)
#--------------------------------
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
##Utils. Fn

get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(base::load(tmp_file))
  tmp_obj<-names(mget(base::load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat<-vroom(paste0(dat_file,cl_res,"/",chromo,".txt"),delim = "\t",col_names = F,trim_ws = T,escape_double = F)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.na(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}

full_bpt_mat<-function(cl_mat,mat_range,res,var){
  
  bin_5kb<-seq(mat_range[1],mat_range[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]
  
  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb
  
  cl_mat$ego_id<-id_conv[as.character(cl_mat$bin.A)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$bin.B)]
  
  #chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=sqrt(-log10(cl_mat$pois.pval)),symmetric = T)
  chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=as.numeric(cl_mat[[var]]))
  return(chr_mat)
}

full_f_mat<-function(cl_mat,res,var){
  
  range_5kb<-range(as.numeric(unique(c(cl_mat$X1,cl_mat$X2))))
  bin_5kb<-seq(range_5kb[1],range_5kb[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]
  
  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb
  
  cl_mat$ego_id<-id_conv[as.character(cl_mat$X1)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$X2)]
  
  #chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=sqrt(-log10(cl_mat$pois.pval)),symmetric = T)
  chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=as.numeric(cl_mat[[var]]),symmetric = T)
  return(chr_mat)
}

#-----------------------------------------
hub_file<-"./data/DAGGER_tbl/trans_res/HMEC_union_trans_res_dagger_tbl.Rda"
top_hub_file<-"./data/DAGGER_tbl/trans_res/HMEC_union_top_trans_res_dagger_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
HiC_dat_folder<-"~/Documents/multires_bhicect/data/HMEC/"

hub_tbl<-get_tbl_in_fn(hub_file) 
top_hub_tbl<-get_tbl_in_fn(top_hub_file) 

tmp_hub_res<-"500kb"
res_hub_candidate<-top_hub_tbl %>% 
  mutate(res=str_split_fixed(node,"_",2)[,1]) %>% 
  filter(res==tmp_hub_res)
hub_res_children_n_tbl<-do.call(bind_rows,map(unique(res_hub_candidate$chr),function(chromo){
  message(chromo)
  tmp_tbl<-res_hub_candidate %>% 
    filter(chr==chromo)
  
  chr_spec_res<-get_tbl_in_fn(paste0(spec_res_file,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  
  chr_hubs_vec<-hub_tbl %>% 
    filter(chr==chromo) %>% 
    dplyr::select(node) %>% 
    unlist
  
  do.call(bind_rows,lapply(tmp_tbl$node,function(candidate_hub){
    candidate_hub_children<-names(which(map_lgl(node_ancestor[chr_hubs_vec],function(i){
      candidate_hub %in% i
    })))
    candidate_hub_children_res<-unique(str_split_fixed(candidate_hub_children,"_",2)[,1])
    
    res_top_children<-vector('list',length(candidate_hub_children_res))
    names(res_top_children)<-candidate_hub_children_res
    for(tmp_res in candidate_hub_children_res){
      
      res_children<-grep(paste0("^",tmp_res),candidate_hub_children,value = T)
      res_top_children[[tmp_res]]<-length(res_children[which(unlist(lapply(res_children,function(x){
        !(any(unlist(lapply(node_ancestor[res_children],function(i){
          x %in% i
        }))))
      })))])
    }
    return(tibble(chr=chromo,hub=candidate_hub,res=candidate_hub_children_res,children=unlist(res_top_children)))
  }))
  
}))
hub_res_children_n_tbl %>%
  mutate(nbin=as.numeric(str_split_fixed(hub,"_",3)[,2])) %>% 
  filter(nbin<=3 & res == "5kb") %>% 
  arrange(desc(children))

#candidate= chr20 500kb_3_3_49000000_50000000  
tmp_hub_chr<-"chr20"
candidate_hub<-"500kb_3_3_49000000_50000000"
chr_spec_res<-get_tbl_in_fn(paste0(spec_res_file,tmp_hub_chr,"_spec_res.Rda"))
chr_bpt<-FromListSimple(chr_spec_res$part_tree)
node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
node_ancestor<-lapply(node_ancestor,'[',-1)

chr_hubs_vec<-hub_tbl %>% 
  filter(chr==tmp_hub_chr) %>% 
  dplyr::select(node) %>% 
  unlist

candidate_hub_children<-names(which(map_lgl(node_ancestor[chr_hubs_vec],function(i){
  candidate_hub %in% i
})))

candidate_hub_children_res<-unique(str_split_fixed(candidate_hub_children,"_",2)[,1])
## Extract top hubs of children at every resolution
res_top_children<-vector('list',length(candidate_hub_children_res))
names(res_top_children)<-candidate_hub_children_res
for(tmp_res in candidate_hub_children_res){
  
  res_children<-grep(paste0("^",tmp_res),candidate_hub_children,value = T)
  res_top_children[[tmp_res]]<-res_children[which(unlist(lapply(res_children,function(x){
    !(any(unlist(lapply(node_ancestor[res_children],function(i){
      x %in% i
    }))))
  })))]
}

mat_range<-range(as.numeric(unique(unlist(chr_spec_res$cl_member[unlist(res_top_children)]))))
res_top_children
# cluster Interaction map by children resolution
tmp_res<-"50kb"
tmp_res_cl<-res_top_children[[tmp_res]]
mat_res_range<-unlist(lapply(mat_range,function(x){
  x - (x %% res_num[tmp_res])
}))

cl_inter_tbl<-do.call(bind_rows,map(seq_along(tmp_res_cl),function(x){
  tmp_cl<-tmp_res_cl[x]
  expand_grid(bin.A=as.numeric(chr_spec_res$cl_member[[tmp_cl]]),bin.B=as.numeric(chr_spec_res$cl_member[[tmp_cl]]),value=x)
}))
cl_mat<-full_bpt_mat(cl_inter_tbl,mat_res_range,res_num[tmp_res],"value")

image(log10(as.matrix(cl_mat)),col=plasma(100))
png(paste0('~/Documents/multires_bhicect/Poster/img/F2/',tmp_hub_chr,"_",candidate_hub,"_",tmp_res,"_cl",'.png'), width =40,height = 40,units = 'mm',type='cairo',res=5000)
par(mar = c(0, 0, 0,0))
plot.new()
image(log10(as.matrix(cl_mat)),col=plasma(100))
dev.off()

#------------------------------
chr_dat<-hic_dat_in(HiC_dat_folder,tmp_res,tmp_hub_chr)
cl_dat<-chr_dat %>% 
  filter(X1 >= mat_res_range[1] & X1 <= mat_res_range[2] & X2 >= mat_res_range[1] & X2 <= mat_res_range[2])
cl_mat<-full_f_mat(cl_dat,res_num[tmp_res],"X3")
image(log10(as.matrix(cl_mat)),col=viridis(100))

png(paste0('~/Documents/multires_bhicect/Poster/img/F2/',tmp_hub_chr,"_",candidate_hub,"_",tmp_res,"_mat",'.png'), width =40,height = 40,units = 'mm',type='cairo',res=5000)
par(mar = c(0, 0, 0,0))
plot.new()
image(log10(as.matrix(cl_mat)),col=viridis(100))
dev.off()
#---------------------------------

tmp_res<-tmp_hub_res
mat_res_range<-unlist(lapply(mat_range,function(x){
  x - (x %% res_num[tmp_res])
}))

cl_inter_tbl<-expand_grid(bin.A=as.numeric(chr_spec_res$cl_member[[candidate_hub]]),bin.B=as.numeric(chr_spec_res$cl_member[[candidate_hub]]),value=2)
cl_mat<-full_bpt_mat(cl_inter_tbl,mat_res_range,res_num[tmp_res],"value")

image(log10(as.matrix(cl_mat)),col=plasma(100))
png(paste0('~/Documents/multires_bhicect/Poster/img/F2/',tmp_hub_chr,"_",candidate_hub,"_",tmp_res,"_cl",'.png'), width =40,height = 40,units = 'mm',type='cairo',res=5000)
par(mar = c(0, 0, 0,0))
plot.new()
image(log10(as.matrix(cl_mat)),col=plasma(100))
dev.off()

#------------------------------
chr_dat<-hic_dat_in(HiC_dat_folder,tmp_res,tmp_hub_chr)
cl_dat<-chr_dat %>% 
  filter(X1 >= mat_res_range[1] & X1 <= mat_res_range[2] & X2 >= mat_res_range[1] & X2 <= mat_res_range[2])
cl_mat<-full_f_mat(cl_dat,res_num[tmp_res],"X3")
image(log10(as.matrix(cl_mat)),col=viridis(100))

png(paste0('~/Documents/multires_bhicect/Poster/img/F2/',tmp_hub_chr,"_",candidate_hub,"_",tmp_res,"_mat",'.png'), width =40,height = 40,units = 'mm',type='cairo',res=5000)
par(mar = c(0, 0, 0,0))
plot.new()
image(log10(as.matrix(cl_mat)),col=viridis(100))
dev.off()
