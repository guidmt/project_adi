create_matrix_meth<-function(work_dir,genes_to_use,nThread=70){
  
  setwd(work_dir)
  
  methylation_directories<-grep(grep(grep(dir(),pattern=".pdf",invert=T,value=T),pattern=".txt",invert=T,value=T),pattern=".png",invert=T,value=T)
  
  # Get a sample to obtain the gene-symbol and probe to use for the analysis
  print("I am selecting the gene-symbol and probe data from a dataset of example")
  setwd(paste(work_dir,methylation_directories[1],sep="/"))
  file_to_read<-dir()
  tab_genes_islands<-fread(file_to_read,select=c(3:6,9,10,11),data.table=T,stringsAsFactors=FALSE,sep="\t")
  tab_genes_islands$Gene_Symbol<-str_split_fixed(string=tab_genes_islands$Gene_Symbol,pattern=";",n=2)[,1]
  tab_genes_islands<-data.frame(tab_genes_islands[Gene_Symbol%in%genes_to_use,])
  
  setwd(work_dir)
  print("start analysis")
  
  registerDoParallel(cores=nThread)
  
  matrixMeth<-foreach(i=1:length(methylation_directories),.combine="cbind",.export =c("fread"),.inorder=T) %dopar% {
    
    print(i)
    
    meth_dir<-methylation_directories[i]
    
    print(meth_dir)
    
    setwd(paste(work_dir,meth_dir,sep="/"))
    
    file_to_read<-dir()
    
    current_sample<-sapply(strsplit(file_to_read,split="\\."),"[[",6)
    #all_samples<-c(all_samples,current_sample)
    #tab_meth<-fread(file_to_read,select=c(3,4,5,6,2,11),data.table=F,nThread=nThread)
    tab_meth<-fread(file_to_read,select=c(6,2),data.table=T,stringsAsFactors=FALSE,sep="\t")
    tab_meth$Gene_Symbol<-str_split_fixed(string=tab_meth$Gene_Symbol,pattern=";",n=2)[,1]
    colnames(tab_meth)[2]<-current_sample
    
    tab_meth2<-data.frame(tab_meth[Gene_Symbol%in%genes_to_use,2])
    
  }
  
  #
  # now add a column with the genes 
  #
  
  matrixMeth2<-cbind(genes=tab_genes_islands,matrixMeth)
  
  return(matrixMeth2)
  
}


