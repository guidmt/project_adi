compareDegsHs_MM<-function(matrix_hs,annotation_for_comp_degs,degs_mouse=ortholog_table_res,lfc=1,pvalue=0.05,output_string){
  
  classes<-annotation_for_comp_degs$status2
  
  PVAL_ALL_GENES<-NULL
  FC_ALL_GENES<-NULL
  LOG2FC_ALL_GENES<-NULL
  MU_num<-NULL
  MU_den<-NULL
  
  for(i in 1:nrow(matrix_hs)){
  
    pvalue_current_gene<-wilcox.test(matrix_hs[i,],classes=classes,correct=FALSE)$p.value
    
    mean_for_groups<-tapply(matrix_hs[i,], classes, mean)
    fc_current_gene<-mean_for_groups["1"]/mean_for_groups["0"]
    log2fc_current_gene<-log2(mean_for_groups["1"]/mean_for_groups["0"])
    
    PVAL_ALL_GENES<-c(PVAL_ALL_GENES,pvalue_current_gene)
    FC_ALL_GENES<-c(FC_ALL_GENES,fc_current_gene)
    LOG2FC_ALL_GENES<-c(LOG2FC_ALL_GENES,log2fc_current_gene)
    
    MU_num<-as.numeric(c(MU_num,mean_for_groups["1"]))
    MU_den<-as.numeric(c(MU_den,mean_for_groups["0"]))
  }
  
  statsDF<-data.frame(genes=rownames(matrix_hs),PVAL_ALL_GENES,MU_num,MU_den,FC_ALL_GENES,LOG2FC_ALL_GENES)
  if(length(which(is.na(statsDF$PVAL_ALL_GENES)))!=0){
  statsDF<-statsDF[-which(is.na(statsDF$PVAL_ALL_GENES)),]
  }else{
  statsDF<-statsDF
  }
  statsDF$padj<-p.adjust(statsDF$PVAL_ALL_GENES,"BH")
  
  # add a label with the status
  statsDF$status<-rep(NA,nrow(statsDF))
  statsDF[statsDF$LOG2FC_ALL_GENES >=lfc & statsDF$PVAL_ALL_GENES<=pvalue,"status"]<-"up"
  statsDF[statsDF$LOG2FC_ALL_GENES <= -lfc & statsDF$PVAL_ALL_GENES<=pvalue,"status"]<-"down"
  statsDF[is.na(statsDF$status),"status"]<-"not_significant"
  
  #
  # Now compare the directions of the genes get in mouse
  #
  
  tcga_up<-as.character(statsDF[statsDF$status=="up",1])
  tcga_down<-as.character(statsDF[statsDF$status=="down",1])
   
  ort_up<-degs_mouse[degs_mouse$status=="up",1]
  ort_up<-as.character(ort_up[!is.na(ort_up)])
  
  ort_down<-degs_mouse[degs_mouse$status=="down",1]
  ort_down<-as.character(ort_down[!is.na(ort_down)])
  
  library(VennDetail)
  ven <- venndetail(list(tcga_up = tcga_up, 
                         tcga_down = tcga_down,
                         ort_up = ort_up,
                         ort_down = ort_down))
  
  venDF<-result(ven, wide = TRUE)
         
  pdf(paste("VennUpSet_DEGS_TCGA_DEGS_Adi",output_string,"pdf",sep="."))
  p1<-plot(ven)
  print(p1)
  p2<-plot(ven, type = "upset")
  print(p2)
  dev.off()
  
  output_xlsx<-paste("VennUpSet_DEGS_TCGA_DEGS_Adi",output_string,"xlsx",sep=".")
  
  write.xlsx(statsDF, file = output_xlsx,sheetName = "TCGA_degs",append = F)
  write.xlsx(venDF, file = output_xlsx,sheetName = "Venn_TCGAdegs_AdiDegs",append = T)
  
}

