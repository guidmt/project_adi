clean_methylation_data<-function(mat_meth,genes_to_use_meth){
  
  ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  tss <- getBM(attributes = c("transcription_start_site", "chromosome_name","transcript_start", "transcript_end", "strand",  "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"),filters = "external_gene_name", values = genes_to_use_meth,mart = ensembl)
  
  mat_meth2<-mat_meth[-which(mat_meth[,7]%in%"."),]
  mat_meth2<-mat_meth2[-which(mat_meth[,7]%in%"."),]
  mat_meth3<-mat_meth2[-which(is.na(mat_meth2[,8])),]
  
  chromosome_cgi<-sapply(strsplit(gsub(mat_meth3[,c(6)],30,pattern="CGI:",replacement=""),split=":"),"[[",1)
  start_cgi<-sapply(strsplit(sapply(strsplit(gsub(mat_meth3[,c(6)],30,pattern="CGI:",replacement=""),split=":"),"[[",2),split="-"),"[[",1)
  end_cgi<-sapply(strsplit(sapply(strsplit(gsub(mat_meth3[,c(6)],30,pattern="CGI:",replacement=""),split=":"),"[[",2),split="-"),"[[",2)
  cgi_df<-data.frame(chromosome_cgi,start_cgi=as.integer(start_cgi),end_cgi=as.integer(end_cgi),CGI=mat_meth3[,c(6)],gene_symbol=mat_meth3[,c(4)],stringsAsFactors = F)
  cgi_df_with_TSS<-merge(tss,cgi_df,by.x="external_gene_name",by.y="gene_symbol")
  cgi_df_with_TSS$overlapped_TSS<-rep(0,nrow(cgi_df_with_TSS))
  cgi_df_with_TSS$putative_genebody<-rep(0,nrow(cgi_df_with_TSS))
  
  #Search for regions that overlaps the TSS: 1) START CGI < TSS & END CGI > TSS + STRAND
  cgi_df_with_TSS$overlapped_TSS[which(cgi_df_with_TSS[,2] > cgi_df_with_TSS$start_cgi & cgi_df_with_TSS[,2] < cgi_df_with_TSS$end_cgi)]<-1
  
  #Search for regions that overlaps the genebody: 1) START TSS < CGI START & END CGI <, strand independent
  idx_genebody1<-which(as.numeric(cgi_df_with_TSS[,2]) < as.numeric(cgi_df_with_TSS$start_cgi))
  idx_genebody2<-which(as.numeric(cgi_df_with_TSS[,2]) < as.numeric(cgi_df_with_TSS$end_cgi))
  gene_body_idx<-intersect(idx_genebody1,idx_genebody2)
  cgi_df_with_TSS$putative_genebody[gene_body_idx]<-1
  
  cgi_df_with_TSS<-cgi_df_with_TSS[-which(cgi_df_with_TSS$putative_genebody %in% 1), ]
  res_cgi_to_use<-unique(cgi_df_with_TSS[which(cgi_df_with_TSS$overlapped_TSS %in% 1), "CGI"])
  mat_meth4<-mat_meth3[mat_meth3[,6]%in%res_cgi_to_use,]
  
  return(mat_meth4)
  
}

combineMethMatrices<-function(lgg_matrix=lgg_matrix,gbm_matrix=gbm_matrix){
  
  if(nrow(lgg_matrix)>nrow(gbm_matrix)){
    
    diff_to_correct<-nrow(lgg_matrix)-nrow(gbm_matrix)
    genes_gbm_meth_df<-data.frame(table(gbm_matrix[,4]))
    genes_lgg_meth_df<-data.frame(table(lgg_matrix[,4]))
    genes_gbm_lgg_all<-cbind(genes_gbm_meth_df,genes_lgg_meth_df)
    genes_to_correct<-genes_gbm_lgg_all[genes_gbm_lgg_all[,2]!=genes_gbm_lgg_all[,4],3]
    
    idx_genes_to_remove_randomly<-sample(which(lgg_matrix[,4]%in%genes_to_correct),1)
    lgg_matrix<-lgg_matrix[-idx_genes_to_remove_randomly,]
    
  }
  
  mat_meth<-cbind(gbm_matrix[,-c(c(1:3),c(5:7))],lgg_matrix[,-c(1:7)])
  colnames(mat_meth)<-gsub(unlist(lapply(strsplit(colnames(mat_meth),split="\\."),FUN=function(x){paste(x[1:3],collapse="-")})),pattern="\\.",replacement="-")
  return(mat_meth) 
}

# mat_meth=combined_meth
# ann_tot=ann_tot
# col_to_use_freq=ann_freq
# genes_to_consider=genes_to_consider
# thr=0.7

compute_meth_levels_for_patient<-function(mat_meth,ann_meth,col_to_use_freq,genes_to_consider,thr){
  
  ann_meth2<-ann_meth[,which(colnames(ann_meth)%in%c("Samples",col_to_use_freq))]
  
  out <- split( ann_meth2 , f = ann_meth2[,colnames(ann_meth2)%in%col_to_use_freq])
  
  df_res_perc<-data.frame()
  
  for(oi in 1:length(out)){

  name_group<-names(out)[oi]
    
  samples_to_use<-out[[oi]][,1]
  mat_meth_sub<-mat_meth[mat_meth[,1]%in%genes_to_consider,colnames(mat_meth)%in%samples_to_use]  
  
  if(length(mat_meth_sub)!=0){
  #the denominator is given by all genes that there are (hyper and hypo in mouse)
  compute_perc<-function(X){length(X[which(X>=thr)])/length(X)}
  
  if(class(mat_meth_sub)=="data.frame"){
  perc_samples_methylated<-apply(mat_meth_sub,2,compute_perc)
  }else{
  perc_samples_methylated<-compute_perc(mat_meth_sub)
  }
  temp_df<-data.frame(name_group,perc_samples_methylated)
  colnames(temp_df)<-c("group","freq")
  
  df_res_perc<-rbind(df_res_perc,temp_df)
  }
  
  }
  
  
  return(df_res_perc)
  
}

compute_statistic<-function(xstat){
  
  stat2<-unlist(lapply(split(xstat,f=xstat$group2),FUN=function(X){wilcox.test(freq~status,data=X)$p.value}))
  return(stat2)
  
}


plot_regions<-function(combined_meth,ann_tot,genes_to_consider_up,genes_to_consider_down,col_to_use,thr,output_pdf="GBM_LGG_hyper_and_hypo_boxplot_on_tcga_filt_wilcox_0.05.pdf"){

    #
    # Hyper regions
    #
  
    hyper_regions<-data.frame()
    
    for(ann_freq in col_to_use){
      print(ann_freq)
      temp_df2<-compute_meth_levels_for_patient(combined_meth,ann_tot,ann_freq,genes_to_consider_up,thr)
      hyper_regions<-rbind(hyper_regions,data.frame(ann_freq,"hyper",temp_df2))
    }

    colnames(hyper_regions)<-c("group1","status","group2","freq")

    #
    # Hypo regions
    #
    
    hypo_regions<-data.frame()
    
    for(ann_freq in col_to_use){
      
      temp_df2<-compute_meth_levels_for_patient(combined_meth,ann_tot,ann_freq,genes_to_consider_down,thr)
      hypo_regions<-rbind(hypo_regions,data.frame(ann_freq,"hypo",temp_df2))
      }
    
    colnames(hypo_regions)<-c("group1","status","group2","freq")
    
    all_hyper_hypo<-rbind(hyper_regions,hypo_regions)
    
    library(forcats)
    
    pdf(output_pdf)
    
    stats_to_use<-round(compute_statistic(all_hyper_hypo),3)
    stats_to_use2<-stats_to_use[stats_to_use<=0.05]
    stats_to_use3<- cbind(rownames(melt(stats_to_use2)),melt(stats_to_use2))
    colnames(stats_to_use3)[1:2]<-c("group_stat","value_stat")
    
    all_hyper_hypo2<-all_hyper_hypo[all_hyper_hypo$group2%in%names(stats_to_use2),]
    all_hyper_hypo3<-merge(all_hyper_hypo2,stats_to_use3,by.x="group2",by.y="group_stat")
    
    p1<-ggplot(all_hyper_hypo3, aes(x = reorder(group2, freq, FUN = mean), y = freq,fill=status)) + geom_boxplot()+ylab("% of methylated regions in a patient")+ylim(0.1,0.5)+facet_grid(group1~., scales = "free", space = "free")+coord_flip()
    print(p1)
    
    grid.newpage()
    
    ss <- grid.table(stats_to_use3)
    print(ss)
    
    dev.off()
    
}

findHumanDMRs<-function(combined_meth,samples_with_mutations,samples_without_mutations){
  
              all_pvalue<-NULL
              all_delta<-NULL
              
              for(i in 1:nrow(combined_meth)){
                
                print(i)
                
                mut<-as.numeric(combined_meth[i,colnames(combined_meth)%in%samples_with_mutations])
                mut_mu<-mean(as.numeric(combined_meth[i,colnames(combined_meth)%in%samples_with_mutations]))
                
                nomut<-as.numeric(combined_meth[i,colnames(combined_meth)%in%samples_without_mutations])
                nomut_mu<-mean(as.numeric(combined_meth[i,colnames(combined_meth)%in%samples_without_mutations]))
                
                all_delta<-c(all_delta,mut_mu-nomut_mu)  
                
                pvalue<-wilcox.test(mut,nomut)$p.value
                all_pvalue<-c(all_pvalue,pvalue)
                
              }
              
              df_stats<-data.frame(combined_meth[,1],all_pvalue,all_delta,p.adjust(all_pvalue,"BH"))
              colnames(df_stats)<-c("genes","all_pvalue","all_delta","padjust")

return(df_stats)

}

correlationMouseWithTCGA<-function(matrix_mouse=res_meth3,matrix_human){
                          
                          colnames(matrix_mouse)[1]<-"gene_symbol"
                          colnames(matrix_human)[1]<-"gene_symbol"
                          
                          mouse_matrix_dt<-setDT(matrix_mouse[!which(is.na(matrix_mouse[,1])),])
                          human_matrix_dt<-setDT(matrix_human)
                          
                          mouse_matrix_rid <- mouse_matrix_dt[, lapply(.SD, mean), by = gene_symbol]
                          human_matrix_rid <- human_matrix_dt[, lapply(.SD, mean), by = gene_symbol]
                          
                          mouse_matrix_rid2<-data.frame(cbind(mouse_matrix_rid[,1],t(apply(mouse_matrix_rid[,-1],1,scale))))
                          colnames(mouse_matrix_rid2)<-colnames(mouse_matrix_rid)
                          
                          human_matrix_rid2<-data.frame(cbind(human_matrix_rid[,1],t(apply(human_matrix_rid[,-1],1,scale))))
                          colnames(human_matrix_rid2)<-colnames(human_matrix_rid)
                          
                          mgiToGeneSymbol = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values =data.frame(mouse_matrix_rid2[,1])[,1], mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
                          mgiToGeneSymbol_commonhs<-mgiToGeneSymbol[mgiToGeneSymbol[,2]%in%data.frame(human_matrix_rid2[,1])[,1],]
                          
                          mouse_matrix3<-merge(mgiToGeneSymbol_commonhs,mouse_matrix_rid2,by.x="MGI.symbol",by.y="gene_symbol")
                          human_matrix3<-merge(mgiToGeneSymbol_commonhs,human_matrix_rid2,by.x="HGNC.symbol",by.y="gene_symbol")
                         
                          mouse_matrix3[is.na(mouse_matrix3)]<-0
                          human_matrix3[is.na(human_matrix3)]<-0
                          
                          #Now compute the correlation between mouse and human
                          res_cor<-data.frame()
                          
                          for(mg_sample in 3:ncol(mouse_matrix3)){
                            
                            all_cor_current_mg_sample<-NULL
                            
                            for(hg_sample in 3:ncol(human_matrix3)){
                              
                              all_cor_current_mg_sample<-c(all_cor_current_mg_sample,cor(mouse_matrix3[,mg_sample],human_matrix3[,hg_sample]))
                              
                            }
                            res_cor<-rbind(res_cor,all_cor_current_mg_sample)
                          }
                          
  
                          rownames(res_cor)<-colnames(mouse_matrix3)[3:ncol(mouse_matrix3)]
                          
                          res_cor2<-cbind(samples=rownames(res_cor),res_cor)
                          
                          res_cor3<-melt(res_cor2)
                          
                      
                          p1<-ggplot(res_cor3, aes(x = samples, y = value,fill=samples)) + geom_boxplot()
} 

dmr_in_mouse=unique(dmrMouse_stats)
dmr_in_human_for_sc=dmr_in_human
  
cross_species_DMRs<-function(dmr_in_human_for_sc,dmr_in_mouse,lfc_meth_mm=1,lfc_meth_hs=0.1){
                    
                    dmr_in_human_for_sc_save<-dmr_in_human_for_sc
  
                    dmr_in_human_for_sc<-dmr_in_human_for_sc[,-c(2,4)]
                    # dmr_in_human_for_pvalue<-dmr_in_human_for_sc[,c(1,4)]
                    
                    dmr_in_human_for_sc<-dmr_in_human_for_sc[!is.na(dmr_in_human_for_sc[,2]),]
                    dmr_in_human_for_sc<-aggregate(dmr_in_human_for_sc[,-1],list(dmr_in_human_for_sc$genes),median)
                    colnames(dmr_in_human_for_sc)[1:2]<-c("gene_symbol","all_delta")

                    colnames(dmr_in_human_for_sc)<-paste(colnames(dmr_in_human_for_sc),"_hs",sep="")
                    colnames(dmr_in_human_for_sc)[1]<-"gene_symbol"
                    colnames(dmr_in_mouse)<-paste(colnames(dmr_in_mouse),"_mm",sep="")
                    colnames(dmr_in_mouse)[1]<-"gene_symbol"
                    
                    mgiToGeneSymbol = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values =dmr_in_mouse[,1], mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
                    
                    mgiToGeneSymbol_common_with_human<-mgiToGeneSymbol[mgiToGeneSymbol[,2]%in%dmr_in_human_for_sc[,1],]
                    
                    dmr_in_human2<-merge(mgiToGeneSymbol_common_with_human,dmr_in_human_for_sc[dmr_in_human_for_sc[,1]%in%mgiToGeneSymbol_common_with_human[,2],],by.x="HGNC.symbol",by.y="gene_symbol")
                    
                    input_for_scatterplot<-merge(dmr_in_human2,dmr_in_mouse,by.x="MGI.symbol",by.y="gene_symbol",all.y=F)
                    
                    input_for_scatterplot$status<-rep("no_sig",nrow(input_for_scatterplot))
                    
                    input_for_scatterplot[input_for_scatterplot$log2FoldChange_mm >0 & input_for_scatterplot$all_delta_hs >0,"status"] <- "hyper_mm/hyper_hs"
                    
                    input_for_scatterplot[input_for_scatterplot$log2FoldChange_mm <0 & input_for_scatterplot$all_delta_hs >0,"status"] <- "hypo_mm/hyper_hs"
                    
                    input_for_scatterplot[input_for_scatterplot$log2FoldChange_mm >0 & input_for_scatterplot$all_delta_hs <0,"status"] <- "hyper_mm/hypo_hs"
                    
                    input_for_scatterplot[input_for_scatterplot$log2FoldChange_mm <0 & input_for_scatterplot$all_delta_hs <0,"status"] <- "hypo_mm-ge/hypo_hs"
                   
                    write.table(input_for_scatterplot,paste("INT_GE_ME",paste(mydir2,current_qval,current_lfc,sep="_"),string_out,"txt",sep="."),sep="\t",row.names=F,quote=F)
                   
                    pdf("test.pdf")
                     
                    require(gg.gap)
                     
                    p<-ggplot(input_for_scatterplot, aes(x=log2FoldChange_mm, y=all_delta_hs)) +
                       geom_point(aes(color=status, fill=status,alpha=0.5),size=3)+scale_shape_manual(values=c(20, 20,1))+scale_color_manual(values=c("blue2","red2","green2", "orange2"))+labs(y="Fold-Change meth HS", x = "Fold-Change meth mm") + theme_classic()+geom_vline(xintercept=-c(-lfc_meth_mm,lfc_meth_mm),linetype="longdash",color="black",size=0.2)+geom_hline(yintercept=c(-lfc_meth_hs,lfc_meth_hs),linetype="longdash",color="black",size=0.2)
                     
                     print(p)
                    
                    dev.off()
                     
}

