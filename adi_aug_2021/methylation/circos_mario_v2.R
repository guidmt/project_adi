library(dplyr)
library(circlize)
library(reshape2)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)

setwd("/home/guidantoniomt/immune_score_manuscript/dictionaries_neo4j_immune")

list_file_relations<-grep(dir(),pattern="relation",value=T)

all_relations<-data.frame()

for(i in list_file_relations){
  
  current_relation<-read.csv(i)
  
  from<-sapply(strsplit(i,split="-"),"[[",1)
  
  to<- sapply(strsplit(sapply(strsplit(i,split="-"),"[[",2),split="_"),"[[",1)
  
    dfinteractions<-data.frame(from=from,to=to,current_relation)
    
    all_relations<-rbind(all_relations,dfinteractions)
  
}

all_relations$sample<-sapply(strsplit(as.character(all_relations$sample),split="-"),"[[",3)

setwd("/home/guidantoniomt/immune_score_manuscript/figs")

circos_mario_fun<-function(all_relations,list_measure="count_corrected",list_function_to_use="sum",use_PFI_scoring=FALSE,norm=FALSE){
  
for(measure in list_measure){
  
  for(func in list_function_to_use){
    
    if(use_PFI_scoring==TRUE){
      
      pdf(paste("CircosInteractionsTMC_PFI_all_ImmuneScoresFoldChange.",func,".",measure,".norm=",norm,".SelfNoSelf.pdf",sep=""))
      
    }else{
      
      pdf(paste("CircosInteractionsTMC_ImmuneScoresFoldChange.",func,".",measure,".norm=",norm,".SelfNoSelf.pdf",sep=""))
      
    }
    
    list_res<-vector(mode="list",2)
    
    labels_to_use<-unique(all_relations$label)
    
    for(ig in  1:length(labels_to_use)){
      
      string_to_use<-labels_to_use[ig]
      
      if(use_PFI_scoring==TRUE){
        
        tab_immune<-read.delim(file="/home/guidantoniomt/immune_score_manuscript/figs/DFI_COAD_ALL_IMMUNE_SCORE_aut_thr.txt")
        tab_immune<-tab_immune[,c("ID","aut_thr")]
        colnames(tab_immune)[2]<-"label"
        tab_immune$label<-ifelse(tab_immune$label=="High_Immune",1,0)
        all_relations2<-merge(tab_immune,all_relations[,-which(colnames(all_relations)%in%("label"))],by.x="ID",by.y="sample")
        
        all_relations2_boxplot<-all_relations2
        all_relations2_boxplot$interaction<-paste(all_relations2_boxplot$from,all_relations2_boxplot$to,sep="_")
        all_relations2_boxplot$label<-ifelse(all_relations2_boxplot$label==1,"High_Immune","Low_Immune")
        
        p<-ggplot(all_relations2_boxplot, aes_string(x="label", y=measure, fill="label"))+
          geom_violin(trim=FALSE)+scale_fill_manual(values=c("red2", "dodgerblue2"))+
          geom_boxplot(width=0.1, fill="white")+ stat_compare_means()+facet_wrap(~interaction, scales = "free",nrow=5)+ylab("Count corrected")+theme_classic()
        
        print(p)
        
      } else {
        
        all_relations2<-all_relations  
        
        all_relations2_boxplot<-all_relations2
        all_relations2_boxplot$interaction<-paste(all_relations2_boxplot$from,all_relations2_boxplot$to,sep="_")
        all_relations2_boxplot$label<-ifelse(all_relations2_boxplot$label==1,"High_Immune","Low_Immune")
        
        p<-ggplot(all_relations2_boxplot, aes_string(x="label", y=measure, fill="label"))+
          geom_violin(trim=FALSE)+scale_fill_manual(values=c("red2", "dodgerblue2"))+
          geom_boxplot(width=0.1, fill="white")+ stat_compare_means()+facet_wrap(~interaction, scales = "free",nrow=5)+ylab("Count corrected")+theme_classic()
        
        print(p)
        
      }
      
      all_relations2<-unique(all_relations2[all_relations2$label==string_to_use,c("from","to",measure)])
      
      if(func=="mean"){
        
        df_relations_stats<-data.frame(all_relations2 %>%
                                         group_by(from, to) %>% 
                                         summarise_each(funs(mean)))
        
      }
      
      if(func=="quantile"){
        
        custom_quantile<-function(x){quantile(x)[4]}
        
        df_relations_stats<-data.frame(all_relations2 %>%
                                         group_by(from, to) %>% 
                                         summarise_each(funs(custom_quantile)))
        
      }
      
      if(func=="sum" & norm == FALSE){
        
        
        df_relations_stats<-data.frame(all_relations2 %>%
                                         group_by(from, to) %>% 
                                         summarise_each(funs(sum)))
      }
      
      if(norm==TRUE){
        
        if(measure=="count_corrected"){    
    
        if(func=="sum"){
          
          #The count corrected are % so convert in this way 
          all_relations2[,3]<-all_relations2[,3]*100
          
          df_relations_stats<-data.frame(all_relations2 %>%
                                           group_by(from, to) %>% 
                                           summarise_each(funs(sum)))
        }
          
        } else{ print("You can do this only with count_corrected")}
        
      }
      
      
      #
      # I need to make symmetrice the table
      #
      symm_matrix<-data.frame(from=df_relations_stats$to,to=df_relations_stats$from,df_relations_stats[,3])
      colnames(symm_matrix)[3]<-measure
      
      df_relations_stats<-rbind(df_relations_stats,symm_matrix)
      
      df_relations_stats_melt<-reshape(df_relations_stats, idvar = "from", timevar = "to", direction = "wide")
      
      df_relations_stats_melt[is.na(df_relations_stats_melt)]<-0
      rownames(df_relations_stats_melt)<-df_relations_stats_melt[,1]
      
      df_relations_stats_melt2<-df_relations_stats_melt[,-1]
      colnames(df_relations_stats_melt2)<-sapply(strsplit(colnames(df_relations_stats_melt2),split="\\."),"[[",2)
      
      df_relations_stats_melt2<-df_relations_stats_melt2[c("tumor","lymph","stroma","misc"),c("tumor","lymph","stroma","misc")]
      
      list_res[[ig]]<-df_relations_stats_melt2
      
    }
    
    names(list_res)<-labels_to_use
    

    if(norm!=TRUE){
      
    mat<-list_res[["1"]]/list_res[["0"]]
    mat[is.na(mat)]<-0
    
    }else{
      
      mat_group1<-list_res[["1"]]
      mat_group1[lower.tri(mat_group1, diag = FALSE)]<-0
      
      mat_group2<-list_res[["0"]]
      mat_group2[lower.tri(mat_group1, diag = FALSE)]<-0
      
      norm_factor_group1<-sum(mat_group1)
      norm_factor_group2<-sum(mat_group2)
        
      list_res[["1"]]<-mat_group1/norm_factor_group1
      list_res[["0"]]<-mat_group2/norm_factor_group2
              
      mat<-list_res[["1"]]/list_res[["0"]]             
      
      # dfbarplot<-rbind(cbind(group="1",melt(as.matrix(list_res[["1"]]))),cbind(group="0",melt(as.matrix(list_res[["0"]]))))
      # dfbarplot$category<-paste(dfbarplot[,2],dfbarplot[,3],sep="_")
      # 
      # ggplot(dfbarplot, aes(x = group, y = value, fill = category)) + 
      #   geom_bar(position = "fill",stat = "identity") +
      #   scale_y_continuous(labels = scales::percent_format())
      
    }
    
    grid.col = c(lymph = "blue2", misc = "darkmagenta", stroma = "orange2", tumor = "red2")
    
    #
    # All interactions self interaction
    #
    mat_rid<-mat
    mat_rid[lower.tri(mat_rid, diag = FALSE)]<-NA
    mat_rid[is.na(mat_rid)]<-0
    
    col_fun = colorRamp2(c(min(mat_rid),0.5,1,1.5,round(max(mat_rid),2)), c("#0066CC","#0080FF","#C0C0C0","#FF9933","#7F00FF"), transparency = 0.5)
    lgd_links = Legend(at = c(min(mat_rid),1,round(max(mat_rid),2)), col_fun = col_fun, title_position = "topleft", title = "Links")
    
    chordDiagram(as.matrix(mat_rid), col = col_fun,grid.col=grid.col,link.border = "black",)
    title(main = "Fold change of TMC HI/LI, all interactions")
    print(draw(lgd_links, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom")))
    
    mat_borders<-mat_rid
    #highlight some interactions
    mat_borders[mat_borders<1.5]<-"NA"
    mat_borders2<-melt(as.matrix(mat_borders))
    if(length(grep(mat_borders2[,3],pattern="NA"))!=0){
    mat_borders2<- mat_borders2[-which(mat_borders2[,3]=="NA"),]
    mat_borders2[,3]<-1
    }

    lty_df = data.frame(mat_borders2[,-3], rep(2,nrow(mat_borders2)))
    lty_df[-which(lty_df[,1]==lty_df[,2]),3]<-1
    lwd_df = data.frame(mat_borders2[,-3],rep(2,nrow(mat_borders2)))
    
    chordDiagram(as.matrix(mat_rid), col = col_fun,grid.col=grid.col,"black",link.lty=lty_df,link.lwd =lwd_df,link.border=mat_borders2)
    title(main = "Fold change of TMC HI/LI, all interactions")
    print(draw(lgd_links, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom")))
    
    #
    # Without self interaction
    #
    
    mat_noself<-mat
    mat_noself[lower.tri(mat_noself, diag = TRUE)]<-NA
    mat_noself[is.na(mat_noself)]<-0
    
    col_fun = colorRamp2(c(min(mat_noself),0.5,1,1.5,round(max(mat_noself),2)), c("#0066CC","#0080FF","#C0C0C0","#FF9933","#7F00FF"), transparency = 0.5)
    lgd_links = Legend(at = c(min(mat_noself),1,round(max(mat_noself),2)), col_fun = col_fun, title_position = "topleft", title = "Links")
    
    chordDiagram(as.matrix(mat_noself), col = col_fun,grid.col=grid.col,link.border = "black")
    title(main = "Fold change of TMC HI/LI, noself interactions")
    print(draw(lgd_links, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom")))
    
    #
    # Self Interaction only
    #
    mat_self<-mat
    mat_self[lower.tri(mat_self, diag = FALSE)]<-NA
    mat_self[upper.tri(mat_self, diag = FALSE)]<-NA
    mat_self[is.na(mat_self)]<-0
    
    col_fun = colorRamp2(c(min(mat_self),0.5,1,1.5,round(max(mat_self),2)), c("#0066CC","#0080FF","#C0C0C0","#FF9933","#7F00FF"), transparency = 0.5)
    lgd_links = Legend(at = c(min(mat_self),1,round(max(mat_self),2)), col_fun = col_fun, title_position = "topleft", title = "Links")
    
    chordDiagram(as.matrix(mat_self), col = col_fun,grid.col=grid.col,link.border = "black")
    title(main = "Fold change of TMC HI/LI, self interactions")
    print(draw(lgd_links, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom")))
    
    #
    # Quality Controls
    #
    
    tme_raw_group1<-melt(cbind(rownames(list_res[["1"]]),data.frame(list_res[["1"]])))
    colnames(tme_raw_group1)<-c("from","to","value")
    
    pc1<-ggplot(tme_raw_group1, aes(to, from)) +
      geom_tile(aes(fill = value)) +
      geom_text(aes(label = round(value, 2))) +
      scale_fill_gradient(low = "white", high = "red")+ggtitle("High Immune: Number of Int.")
    
    tme_raw_group2<-melt(cbind(rownames(list_res[["0"]]),data.frame(list_res[["0"]])))
    colnames(tme_raw_group2)<-c("from","to","value")
    
    pc2<-ggplot(tme_raw_group2, aes(to, from)) +
      geom_tile(aes(fill = value)) +
      geom_text(aes(label = round(value, 2))) +
      scale_fill_gradient(low = "white", high = "red")+ggtitle("Low Immune: Number of Int.")
    
    fold_change_tcm<-melt(cbind(rownames(mat),data.frame(mat)))
    colnames(fold_change_tcm)<-c("from","to")
    
    pc3<-ggplot(melt(fold_change_tcm), aes(to, from)) +
      geom_tile(aes(fill = value)) +
      geom_text(aes(label = round(value, 2))) +
      scale_fill_gradient(low = "white", high = "blue")+ggtitle("Number of Int.")
    
    print(plot_grid(pc1,pc2,pc3))
    
    dev.off()
    
  }
  
}

}

circos_mario_fun(all_relations,list_measure="count_corrected",list_function_to_use="sum",use_PFI_scoring=FALSE,norm=TRUE)
circos_mario_fun(all_relations,list_measure="count_corrected",list_function_to_use="sum",use_PFI_scoring=FALSE,norm=FALSE)

circos_mario_fun(all_relations,list_measure="count_corrected",list_function_to_use="sum",use_PFI_scoring=TRUE,norm=TRUE)
circos_mario_fun(all_relations,list_measure="count_corrected",list_function_to_use="sum",use_PFI_scoring=TRUE,norm=FALSE)


