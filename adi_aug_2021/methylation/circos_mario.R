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
  
  if(from!=to){
    
  dfinteractions<-data.frame(from=from,to=to,current_relation)

  all_relations<-rbind(all_relations,dfinteractions)
  
  }

}

all_relations$sample<-sapply(strsplit(as.character(all_relations$sample),split="-"),"[[",3)

setwd("/home/guidantoniomt/immune_score_manuscript/figs")

list_measure<-c("count_corrected","count")

for(measure in list_measure){

list_function_to_use<-c("mean","sum")

use_PFI_scoring<-TRUE
  
for(func in list_function_to_use){

if(use_PFI_scoring==TRUE){
  
pdf(paste("CircosInteractionsTMC_PFI_all_ImmuneScoresFoldChange.",func,".",measure,".pdf",sep=""))

}else{
  
pdf(paste("CircosInteractionsTMC_ImmuneScoresFoldChange.",func,".",measure,".pdf",sep=""))
  
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
      geom_boxplot(width=0.1, fill="white")+ stat_compare_means()+facet_wrap(~interaction, scales = "free")+ylab("Count corrected")+theme_classic()
    
    print(p)
    
  } else {
    
  all_relations2<-all_relations  
  
  all_relations2_boxplot<-all_relations2
  all_relations2_boxplot$interaction<-paste(all_relations2_boxplot$from,all_relations2_boxplot$to,sep="_")
  all_relations2_boxplot$label<-ifelse(all_relations2_boxplot$label==1,"High_Immune","Low_Immune")
  
  p<-ggplot(all_relations2_boxplot, aes_string(x="label", y=measure, fill="label"))+
    geom_violin(trim=FALSE)+scale_fill_manual(values=c("red2", "dodgerblue2"))+
    geom_boxplot(width=0.1, fill="white")+ stat_compare_means()+facet_wrap(~interaction, scales = "free")+ylab("Count corrected")+theme_classic()
  
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
  
  if(func=="sum"){
    
    custom_quantile<-function(x){quantile(x)[4]}
    
    df_relations_stats<-data.frame(all_relations2 %>%
                                     group_by(from, to) %>% 
                                     summarise_each(funs(sum)))
    
  }
  
  df_relations_stats_melt<-reshape(df_relations_stats, idvar = "from", timevar = "to", direction = "wide")
  
  df_relations_stats_melt[is.na(df_relations_stats_melt)]<-0
  rownames(df_relations_stats_melt)<-df_relations_stats_melt[,1]
  
  df_relations_stats_melt2<-df_relations_stats_melt[,-1]
  colnames(df_relations_stats_melt2)<-sapply(strsplit(colnames(df_relations_stats_melt2),split="\\."),"[[",2)
  
  # chordDiagram(as.matrix(df_relations_stats_melt2), transparency = 0.5,annotationTrack = c("name", "grid"))
  # title(main = string_to_use)
  
  list_res[[ig]]<-df_relations_stats_melt2
  
}

names(list_res)<-labels_to_use

grid.col = c(lymph = "blue2", misc = "darkmagenta", stroma = "orange2", tumor = "red2")

mat<-list_res[["1"]]/list_res[["0"]]
mat[is.na(mat)]<-0

col_fun = colorRamp2(c(min(mat),0.5,1,1.5,round(max(mat),2)), c("#0066CC","#0080FF","#C0C0C0","#FF9933","#7F00FF"), transparency = 0.5)
lgd_links = Legend(at = c(min(mat),1,round(max(mat),2)), col_fun = col_fun, title_position = "topleft", title = "Links")

chordDiagram(as.matrix(mat), col = col_fun,grid.col=grid.col,link.border = "black")
title(main = "Fold change of TMC interactions between High Immune and Low Immune")
draw(lgd_links, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

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

