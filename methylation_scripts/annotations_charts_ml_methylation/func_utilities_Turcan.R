prepareDataTCGAandTurcan<-function(hs_meth_in_mouse_gbm,hs_meth_in_mouse_lgg,genes_to_use,ann3,rand=F){
  
  datasets_parsing<-vector(mode="list",2)
  
  #
  # Prepare TCGA data
  #
  
  hs_meth_in_mouse_gbm_filt<-clean_methylation_data(mat_meth=hs_meth_in_mouse_gbm,genes_to_use_meth=genes_to_use)
  hs_meth_in_mouse_lgg_filt<-clean_methylation_data(mat_meth=hs_meth_in_mouse_lgg,genes_to_use_meth=genes_to_use)
  
  ann_tot<-ann3
  
  if(rand==F){
  combined_meth<-combineMethMatrices(lgg_matrix=hs_meth_in_mouse_lgg_filt,gbm_matrix=hs_meth_in_mouse_gbm_filt)
  }else{
  combined_meth<-combineMethMatricesRandomAnalysis(lgg_matrix=hs_meth_in_mouse_lgg_filt,gbm_matrix=hs_meth_in_mouse_gbm_filt)
  }
  colnames(combined_meth)[1]<-"HGNC.symbol"
  
  human_tab<-combined_meth
  
  input_for_ml<-t(human_tab[,-1])
  colnames(input_for_ml)<-human_tab[,1]
  input_for_ml2<-data.frame(ID=rownames(input_for_ml),input_for_ml)
  input_for_ml3<-merge(ann_tot,input_for_ml2,by.x="Case",by.y="ID",all.x=F,all.y=F)
  colnames(input_for_ml3)<-gsub(colnames(input_for_ml3),pattern="\\.",replacement="-")
  
  colnames(ann3)[1]<-"Samples"
  
  input_ml4<-input_for_ml3[,-c(2:19)]
  ann3$IDH_TOT<-rep("0",nrow(ann3))
  ann3[which(ann3$IDH1==1 | ann3$IDH2==1),"IDH_TOT"]<-1
  
  input_ml5<-merge(ann3,input_ml4,by.x="Samples",by.y="Case")
  
  #
  # Prepare Turcan data
  #
  
  setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/public_data")
  
  tab_input<-fread("GSE30338_81clinic_53cellline_rawbeta.txt.gz",data.table=F,fill=T)
  
  tab_input2<-tab_input[,colnames(tab_input)%in%c("Gene Symbol",grep(colnames(tab_input),pattern="^s",value=T))]
  tab_input3<-cbind(gene_symbol=tab_input2[,ncol(tab_input2)],tab_input2[,-ncol(tab_input2)])
  
  annotation_from_geo<-read.delim("GSE30338_annotation_from_geo.txt",stringsAsFactors = F)
  
  ann_paper <- data.frame(read_excel("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/public_data/Table 5a.xls"))
  colnames(ann_paper)<-ann_paper[1,]
  ann_paper2<-ann_paper[-1,]
  
  annotation_from_geo$Tumour_id<-as.character(annotation_from_geo$Tumour_id)
  annotation_final<-merge(annotation_from_geo,ann_paper2,by.x="Tumour_id",by.y="ID number")
  
  tab_input3<-tab_input3[-which(tab_input3$gene_symbol==""),]
  
  gs<-strsplit(as.character(tab_input3[,1]),split="\\;")
  gs2<-lapply(gs,FUN=unique)

  gs3<-unlist(lapply(gs2,FUN=function(X){paste0(X,collapse=";")}))
  tab_input3$gene_symbol<-gs3
  
  check_illumina<-functionCheckGeneMet(gs2,list_genes=genes_to_use)
  
  mat_turcan<-data.table(tab_input3[check_illumina[[1]],])
  mat_turcan_agg <- data.frame(mat_turcan[, lapply(.SD, mean), by = gene_symbol])
  
  mat_turcan2<-t(mat_turcan_agg[,-1])
  mat_turcan3<-data.frame(ID=rownames(mat_turcan2),mat_turcan2)
  mat_turcan4<-merge(annotation_final[,c(2,4)],mat_turcan3,by.x="sampleID",by.y="ID")
  mat_turcan4[,2]<-ifelse(mat_turcan4[,2]=="CIMP+","IDH_yes","IDH_no")
  colnames(mat_turcan4)[2]<-"status"
  colnames(mat_turcan4)[3:ncol(mat_turcan4)]<-mat_turcan_agg[,1]
  
  datasets_parsing[[1]]<-input_ml5
  datasets_parsing[[2]]<-mat_turcan4
  
  return(datasets_parsing)
  
}

functionCheckGeneMet<-function(gs2,list_genes=c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2])){
  
  index_illumina<-NULL
  custom_list<-NULL
  
  # index_custom<-data.frame()
  
  # res_match_illumina<-vector(mode="list",2)
  
  for(i in 1:length(gs2)){
    
    print(i)
    
    #i get all the gene symbol annotated with one CGI,
    #if there is only 1 gene string_genes is the current genes
    #see if the gene annotated is present in your gene list
    
    if(length(gs2[[i]])==1){
      
      string_genes<- paste("^",gs2[[i]],"$",sep="")
      
      ovp<-length(grep(list_genes,pattern=string_genes))
      
    } else{
      
      #i get all the gene symbol annotated with one CGI,
      #if there is more than 1 gene string_genes  contains all the genes for the given cgi
      #see if the gene annotated is present in your gene list
      
      string_genes<- paste(paste("^",gs2[[i]],"$",sep=""),collapse="|")
      
      ovp<-length(grep(list_genes,pattern=string_genes))
      
    }
    
    #if there is any overlap save the index of the i-th cgi
    if(ovp>=1){
      
      index_illumina<-c(index_illumina,i)
      indeces_grep<-paste(grep(list_genes,pattern=string_genes),collapse=";")
      custom_list<-c(custom_list,indeces_grep)
    }
    
  }
  
  res_check_match_illumina<-vector(mode="list",2)
  
  res_check_match_illumina[[1]]<-index_illumina
  res_check_match_illumina[[2]]<-custom_list
  
  return(res_check_match_illumina) 
}
