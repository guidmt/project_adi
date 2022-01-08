library(data.table)

# 
#  Script that take the DMRs genes and compute the % of samples with hyper or hypo methylation across different TCGA-subgroups
# 

dmrMouse<-read.delim(file="/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/annotation_dmrs_adi_pval0.05_lfc1.txt",stringsAsFactors = F)

dmrMouse_rid<-unique(dmrMouse[,c(1,2,3,19)])

# 
#  Extract the beta values in mouse 
# 

library(data.table)

setwd("/mnt/data/lab/gmt_data/data_brain/adi/CytosineReports")

list_meth_file<-dir()

list_samples<-vector(mode='list',length(list_meth_file))

for(i in 1:length(list_meth_file)){
  
  print(i)
  
  print("reading data")
  meth_table_current_mouse<-fread(file=list_meth_file[i])
  colnames(meth_table_current_mouse)<-c("chromosome","position","strand","count_methylated","count_unmethylated","c_context","trinucleotide_context")
  
  beta_current_values<-data.frame()
  
  #select only the sites in which the gens were identified DEGs
  for(row_ann in 1:nrow(dmrMouse_rid)){
    
    print(paste("Genes:",row_ann))
    
    chrom<-dmrMouse_rid[row_ann,1]
    start<-dmrMouse_rid[row_ann,2]
    end<-dmrMouse_rid[row_ann,3]
    genes<-dmrMouse_rid[row_ann,4]
    
    temp_chr<-meth_table_current_mouse[chromosome==chrom]
    
    test<-data.frame(temp_chr[ position >=start & position <= end ])
    
    if(nrow(test)>=1){
      
      coordinates<-paste(range(test[,1]),collapse="-")
      # https://www.biostars.org/p/445286/
      beta_values<-sum(test$count_methylated)/(sum(test$count_methylated)+sum(test$count_unmethylated))
      
      df_beta<-data.frame(gene_symbol=genes,beta_values)
      
      beta_current_values<-rbind(beta_current_values,df_beta)
      
    }
    
  }
  
  
  list_samples[[i]]<-beta_current_values
  
}

res_meth<-do.call(cbind,list_samples)

res_meth2<-res_meth[,grep(colnames(res_meth),pattern="beta_values")]
colnames(res_meth2)<-sapply(strsplit(list_meth_file,split="\\."),"[[",1)

res_meth3<-as.data.frame(unique(cbind(genes=as.character(res_meth[,1]),res_meth2)))

remove_inf<-function(x){
  idx<-which(!is.finite(x))
  x[idx]<-0
  return(x)
}

res_meth3[,-1]<-sapply(res_meth3[,-1],remove_inf)

colnames(res_meth3)[1]<-"gene_symbol"

save("Mouse_beta_matrix_adi_pval0.05_lfc1.RData")

# 
#  Compute % of samples with genes that are hyper and hypo-methylated in mouse and the two groups then you can create a boxplot of this
# 

