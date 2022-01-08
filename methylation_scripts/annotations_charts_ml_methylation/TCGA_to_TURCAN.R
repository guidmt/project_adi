library("biomaRt")
library("data.table")
library(plyr)
library(caret)
library(cowplot)
library(patchwork)
library(data.table)
library(minfi)
library("readxl")
library(factoextra)
library(sva)
library(LiblineaR)

# 
# Download all the mouse genes 
# 

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mmGenes<-getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), mart = mouse)

# setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")
# load("Mouse_beta_matrix_adi_pval0.05_lfc1.RData")

source("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/clean_methylation_data.R")

#
# Read the DMRs regions of mouse
#

dmrMouse<-read.delim(file="/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/annotation_dmrs_adi_pval0.05_lfc1.txt",stringsAsFactors = F)
dmrMouse_rid<-unique(dmrMouse[,c(1,2,3,19,17)])
dmrMouse_stats<-dmrMouse[,c(19,13,14,17,18)]

#
# Convert mouse genes in human genes
#

httr::set_config(httr::config(ssl_verifypeer = FALSE))

human = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
mouse = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")

dmrs_up_mm<-unique(dmrMouse[dmrMouse$status=="DMRs-up",19])
dmrs_down_mm<-unique(dmrMouse[dmrMouse$status=="DMRs-down",19])

hs_dmrs_up_from_mm = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dmrs_up_mm , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

hs_dmrs_down_from_mm = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dmrs_down_mm , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

#
# Prepare the TCGA human methylation data and annotation (Illumina InfiniumHumanMethylation450 platform)
#

#look the script: Methylation_barplot_mouse_human.R for pre-processing methylation data
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")
# see script: Methylation_barplot_mouse_human.R to see how I have created the RData

load("GBM_Mouse_beta_matrix_adi_pval0.05_lfc1.RData")
load("LGG_Mouse_beta_matrix_adi_pval0.05_lfc1.RData")


hs_meth_in_mouse_gbm_filt<-clean_methylation_data(mat_meth=hs_meth_in_mouse_gbm,genes_to_use_meth=c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2]))

hs_meth_in_mouse_lgg_filt<-clean_methylation_data(mat_meth=hs_meth_in_mouse_lgg,genes_to_use_meth=c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2]))

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human")
load("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/annotation_gbm_lgg_July2021.RData")

ann_tot<-ann3

#  Step 4: Combine all methylation
combined_meth<-combineMethMatrices(lgg_matrix=hs_meth_in_mouse_lgg_filt,gbm_matrix=hs_meth_in_mouse_gbm_filt)
colnames(combined_meth)[1]<-"HGNC.symbol"
  
# colnames(hs_meth_in_mouse_lgg_filt)<-gsub(unlist(lapply(strsplit(colnames(hs_meth_in_mouse_lgg_filt),split="\\."),FUN=function(x){paste(x[1:3],collapse="-")})),pattern="\\.",replacement="-")
# colnames(hs_meth_in_mouse_gbm_filt)<-gsub(unlist(lapply(strsplit(colnames(hs_meth_in_mouse_gbm_filt),split="\\."),FUN=function(x){paste(x[1:3],collapse="-")})),pattern="\\.",replacement="-")


#
# Prepare the input for machine learning
#

require(caret)

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")
source("libLinear_mouse_hsApril.R")

# meth_data_mouse<-setDT(res_meth3)
# meth_data_mouse2 <- data.frame(meth_data_mouse[, lapply(.SD, mean), by = gene_symbol])
# 
# httr::set_config(httr::config(ssl_verifypeer = FALSE))
# 
# human = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
# mouse = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")
# 
# #convert the mouse genes in human genes
# genes_mm_map_hs = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = meth_data_mouse2[,1] , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
# 
# #Measure the average Beta values for the genes in homo sapiens
# colnames(combined_meth)[1]<-"gene_symbol"
# meth_data_hs<-setDT(combined_meth)
# meth_data_hs2 <- data.frame(meth_data_hs[, lapply(.SD, mean), by = gene_symbol])
# 
# #get the common genes between human and mouse
# common_genes_hs_mm<-intersect(meth_data_hs2[,1],genes_mm_map_hs[,2])
# 
# mouse_genes_to_consider<-genes_mm_map_hs[genes_mm_map_hs[,2]%in%common_genes_hs_mm,1]
# human_genes_to_consider<-genes_mm_map_hs[genes_mm_map_hs[,2]%in%common_genes_hs_mm,2]
# 
# #Create the mouse table
# meth_data_mouse3<-merge(genes_mm_map_hs,meth_data_mouse2[meth_data_mouse2[,1]%in%mouse_genes_to_consider,],by.x="MGI.symbol",by.y="gene_symbol")
# #Create the mouse table
# meth_data_hs3<-meth_data_hs2[meth_data_hs2[,1]%in%human_genes_to_consider,]

# Now merge the two data.frame because is the only way to get the same number of columns
# Warning: i will not use more the mouse gene symbol but only the human gene symbol to be consistent

# merged_matrix<-merge(meth_data_mouse3,meth_data_hs3,by.x="HGNC.symbol",by.y="gene_symbol")

# Re-extract the human data and prepare the input for ML
# idx1<-grep(colnames(merged_matrix),pattern="TCGA")
# idx2<-grep(colnames(merged_matrix),pattern="HGNC.symbol")
# human_tab<-merged_matrix[,c(idx2,idx1)]
human_tab<-combined_meth
# colnames(human_tab)<-gsub(colnames(human_tab),pattern="\\.",replacement="-")

input_for_ml<-t(human_tab[,-1])
colnames(input_for_ml)<-human_tab[,1]
input_for_ml2<-data.frame(ID=rownames(input_for_ml),input_for_ml)
input_for_ml3<-merge(ann_tot,input_for_ml2,by.x="Case",by.y="ID",all.x=F,all.y=F)
colnames(input_for_ml3)<-gsub(colnames(input_for_ml3),pattern="\\.",replacement="-")

# idx1<-grep(colnames(merged_matrix),pattern="TCGA",invert=T)
# idx2<-grep(colnames(merged_matrix),pattern="HGNC.symbol")
# mouse_tab<-merged_matrix[,c(idx2,idx1)][,-c(2:4)] #remove columns with gene-symbol and control EC

# input_mm<-t(mouse_tab[,-1])
# colnames(input_mm)<-mouse_tab[,1]
# input_mm2<-data.frame(ID=c("mut","mut","mut","wt","wt","wt","wt"),input_mm)
# colnames(input_mm2)[2:ncol(input_mm2)]<-gsub(colnames(input_mm2)[2:ncol(input_mm2)],pattern="\\.",replacement="-")
# colnames(input_mm2)[1]<-"status_mouse"
# rownames(input_mm2)<-colnames(mouse_tab)[-1]

list_ann<-c("DNA_methylation_clusters","IDH_TOT")
source("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/libLinear_mouse_hsApril.R")

#
# Add the annotation with the input of machine learning
#

load("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/annotation_gbm_lgg_July2021.RData")
colnames(ann3)[1]<-"Samples"

input_ml4<-input_for_ml3[,-c(2:19)]
ann3$IDH_TOT<-rep("0",nrow(ann3))
ann3[which(ann3$IDH1==1 | ann3$IDH2==1),"IDH_TOT"]<-1

#
# input_ml5 is the input of machine learning for TCGA
#

input_ml5<-merge(ann3,input_ml4,by.x="Samples",by.y="Case")

######################################
# Prepare the turcan data (Illumina InfiniumHumanMethylation450 platform)
######################################

# define a function to select the Illumina Probe IDs related with the DMRs regions
functionCheckGeneMet<-function(list_meth=gs2,list_genes=c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2])){
  
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

dmrMouse<-read.delim(file="/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/annotation_dmrs_adi_pval0.05_lfc1.txt",stringsAsFactors = F)

require("biomaRt")
#Warning-use the last version of Ensembl to convert mgi symbol to hg symbol - call useEnsembl to use a different version
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

dmrs_up_mm<-unique(dmrMouse[dmrMouse$status=="DMRs-up",19])
hs_dmrs_up_from_mm = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dmrs_up_mm , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

dmrs_down_mm<-unique(dmrMouse[dmrMouse$status=="DMRs-down",19])
hs_dmrs_down_from_mm = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dmrs_down_mm , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)


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

#
# Exploratory analysis
#

tab_input3<-tab_input3[-which(tab_input3$gene_symbol==""),]

gs<-strsplit(as.character(tab_input3[,1]),split="\\;")
gs2<-lapply(gs,FUN=unique)

gs3<-unlist(lapply(gs2,FUN=function(X){paste0(X,collapse=";")}))
tab_input3$gene_symbol<-gs3

check_illumina<-functionCheckGeneMet(gs2,c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2]))

mat_turcan<-data.table(tab_input3[check_illumina[[1]],])
mat_turcan_agg <- data.frame(mat_turcan[, lapply(.SD, mean), by = gene_symbol])

mat_turcan2<-t(mat_turcan_agg[,-1])
mat_turcan3<-data.frame(ID=rownames(mat_turcan2),mat_turcan2)
mat_turcan4<-merge(annotation_final[,c(2,4)],mat_turcan3,by.x="sampleID",by.y="ID")
mat_turcan4[,2]<-ifelse(mat_turcan4[,2]=="CIMP+","IDH_yes","IDH_no")
colnames(mat_turcan4)[2]<-"status"
colnames(mat_turcan4)[3:ncol(mat_turcan4)]<-mat_turcan_agg[,1]
  
#
# Use the turcan datasets as a validation set for TCGA only in LGG
#

tcga_dataset<-input_ml5[which(input_ml5$type%in%"LGG"),20:ncol(input_ml5)]
tcga_dataset[,1]<-as.factor(tcga_dataset[,1])
colnames(tcga_dataset)[1]<-"status"
tcga_dataset[is.na(tcga_dataset)]<-0
tcga_dataset$status<-ifelse(tcga_dataset$status==1,"IDH_yes","IDH_no")
  
common_genes_datasets<-intersect(colnames(tcga_dataset)[-1],colnames(mat_turcan4)[-1])
tcga_dataset2<-tcga_dataset[,which(colnames(tcga_dataset)%in%c("status",common_genes_datasets))]
mat_turcan5<-mat_turcan4[,which(colnames(mat_turcan4)%in%c("status",common_genes_datasets))]

library(caret)

ctrl <- trainControl(method="cv",
                     number=5,
                     classProbs=T,
                     savePredictions=T,
                     sampling = "down",
                     summaryFunction=twoClassSummary,
                     allowParallel=TRUE)
print("train")
rfmodel <- train(status ~ .,
                 data=tcga_dataset2,
                 method="ranger",
                 trControl=ctrl)

rf_pred <- predict(rfmodel, mat_turcan5[,-1])

overall_per<-confusionMatrix(rf_pred, as.factor(mat_turcan4$status),positive="IDH_yes")$overall
byclass_per<-confusionMatrix(rf_pred, as.factor(mat_turcan4$status),positive="IDH_yes")$byClass

write.table(overall_per,file="TCGA_to_TURCAN_NoCombat_overall_performance.txt",row.names=F,quote=F,sep="\t")
write.table(byclass_per,file="TCGA_to_TURCAN_NoCombat_BinaryClassPerformance.txt",row.names=F,quote=F,sep="\t")

#
# PCA to check the data
#

input_pca<-rbind.fill(tcga_dataset,mat_turcan5)
input_pca[is.na(input_pca)]<-0
res.pca <- prcomp(input_pca[,-1], scale = TRUE)

library(factoextra)

pdf("PCA_Turcan_withTCGA.DMRs.mouse.pdf")
p<-fviz_pca_ind(res.pca,
                col.ind = input_pca[,1], # color by groups
                palette = c("#00AFBB",  "#FC4E07"),
                legend.title = "Groups",
                label="none",
                pointsize = 2.5
)
print(p)
p2<-fviz_pca_ind(res.pca,
                col.ind = c(rep("TCGA",nrow(tcga_dataset)),rep("Turcan",nrow(mat_turcan5))), # color by groups
                palette = c("midnightblue",  "red4"),
                legend.title = "Groups",
                label="none",
                pointsize = 2.5
)
print(p2)
dev.off()

#
# Remove batch - effects
#

all_status<-rbind.fill(tcga_dataset2,mat_turcan5)[,1]
meth_all<-t(as.matrix(rbind.fill(tcga_dataset2,mat_turcan5)[,-1]))
meth_all[is.na(meth_all)]<-0

groups<-as.factor(c(rep("TCGA",nrow(tcga_dataset2[,-1])),rep("Turcan",nrow(mat_turcan5[,-1]))))

combat_meth = ComBat(dat=meth_all, batch=groups)

combat_meth2<-data.frame(Groups=all_status,t(combat_meth))

res.pca_combat <- prcomp(combat_meth2[,-1], scale = TRUE)


pdf("PCA_PostCombat_Turcan_withTCGA.DMRs.mouse.pdf")
p<-fviz_pca_ind(res.pca_combat,
                col.ind = combat_meth2[,1], # color by groups
                palette = c("#00AFBB",  "#FC4E07"),
                legend.title = "Groups",
                label="none",
                pointsize = 2.5
)
print(p)
p2<-fviz_pca_ind(res.pca_combat,
                 col.ind = groups, # color by groups
                 palette = c("midnightblue",  "red4"),
                 legend.title = "Groups",
                 label="none",
                 pointsize = 2.5
)
print(p2)
dev.off()

#
# Use the machine learning to see the predictions between TCGA and Turcan
#

combat_meth2<-data.frame(Groups=groups,IDH_status=all_status,t(combat_meth))

train_tcga_combat<-combat_meth2[combat_meth2$Groups%in%"TCGA",-1]
test_tcga_combat<-combat_meth2[combat_meth2$Groups%in%"Turcan",-1]

ctrl <- trainControl(method="cv",
                     number=5,
                     classProbs=T,
                     savePredictions=T,
                     sampling = "down",
                     summaryFunction=twoClassSummary,
                     allowParallel=TRUE)
print("train")
rfmodel <- train(IDH_status ~ .,
                 data=train_tcga_combat,
                 method="ranger",
                 trControl=ctrl)

rf_pred <- predict(rfmodel, test_tcga_combat[,-1])

overall_per<-confusionMatrix(rf_pred, as.factor(mat_turcan4$status),positive="IDH_yes")$overall
byclass_per<-confusionMatrix(rf_pred, as.factor(mat_turcan4$status),positive="IDH_yes")$byClass

write.table(overall_per,file="TCGA_to_TURCAN_Combat_overall_performance.txt",row.names=F,quote=F,sep="\t")
write.table(byclass_per,file="TCGA_to_TURCAN_Combat_BinaryClassPerformance.txt",row.names=F,quote=F,sep="\t")


# rf_pred <- predict(rfmodel, test_tcga_combat[,-1],type="prob")
# 
# rf.ROC <- roc(predictor=rf_pred$IDH_yes,response=test_tcga_combat$IDH_status)
# confidence_interval<-ci.auc(test_tcga_combat$IDH_status, rf_pred$IDH_yes, conf.level = 0.9)
# ci_string<-paste(round(as.numeric(confidence_interval),3)[2],"(",round(as.numeric(confidence_interval),3)[1],"-",round(as.numeric(confidence_interval),3)[3],")",sep="")
# 
# pAUCCI<-ggroc(rf.ROC) + theme_minimal() +
#   geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + coord_equal() +
#   geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha= 0.2) + ggtitle(ci_string)

#
# Liblinear does not work really properly on these data
#
# co=heuristicC(as.matrix(train_tcga_combat[,-1]))
# model <- LiblineaR(data=train_tcga_combat[,-1],target=train_tcga_combat[,1],cost=co,type=0)
# predictions<-LiblineaR:::predict.LiblineaR(model, test_tcga_combat[,-1],proba=T)$predictions
