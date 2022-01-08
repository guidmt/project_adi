
create_matrix_meth<-function(methylation_human_dir,genes_to_use,nThread=60){
  
  current_dir<-paste("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/GDCdata",methylation_human_dir,"harmonized/DNA_Methylation/Methylation_Beta_Value",sep="/")
  
  setwd(current_dir)
  
  methylation_directories<-grep(grep(grep(dir(),pattern=".pdf",invert=T,value=T),pattern=".txt",invert=T,value=T),pattern=".png",invert=T,value=T)
  
  matrix_meth<-vector(mode="list",length(methylation_directories))
  
  name_samples<-NULL
  
  for(i in 1:length(methylation_directories)){
    
    print(i)
    
    meth_dir<-methylation_directories[i]
    
    setwd(paste(current_dir,meth_dir,sep="/"))
    
    file_to_read<-dir()
    current_sample<-sapply(strsplit(file_to_read,split="\\."),"[[",6)
    name_samples<-c(name_samples,current_sample)
    
    tab_meth<-fread(file_to_read,select=c(3,4,5,6,2,11),data.table=F,nThread=nThread)
    
    tab_meth$Gene_Symbol<-str_split_fixed(string=tab_meth$Gene_Symbol,pattern=";",n=2)[,1]

    tab_meth2<-tab_meth[which(tab_meth$Gene_Symbol %in% genes_to_use),]
    
    matrix_meth[[i]]<-tab_meth2
    
    setwd(current_dir)
    
  }
  
  matrix_meth2<-do.call(cbind,matrix_meth)
  
  idx_column_to_select<-grep(colnames(matrix_meth2),pattern="Beta_value")
  
  matrix_meth3<-matrix_meth2[,idx_column_to_select]
  
  colnames(matrix_meth3)<-name_samples
  
  matrix_meth4<-data.frame(gene_symbol=matrix_meth[[1]][,c(4,6)],matrix_meth3)
  
  matrix_meth5<-matrix_meth4[,grep(sapply(strsplit(colnames(matrix_meth4),split="\\."),"[",4),pattern=paste("10","11","12","13","14",collapse="|"),invert=T)]
  
  return(matrix_meth5)
}
