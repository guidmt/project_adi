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
# setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")
# load("Mouse_beta_matrix_adi_pval0.05_lfc1.RData")

source("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/clean_methylation_data.R")
source("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/func_utilities_Turcan.R")
load("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/annotation_gbm_lgg_July2021.RData")

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

genes_mouse_to_use<-c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2])

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")

load("GBM_Mouse_beta_matrix_adi_pval0.05_lfc1.RANDOMGENES.RData")
load("LGG_Mouse_beta_matrix_adi_pval0.05_lfc1.RANDOMGENES.RData")

size_rand_genes<-c(length(genes_mouse_to_use),round(length(genes_mouse_to_use)/seq(2,20,by=2)))

all_stats<-data.frame()

for(i in size_rand_genes){

genes_mouse_to_use_rand<-sample(setdiff(x=unique(intersect(hs_meth_in_mouse_lgg$genes.Gene_Symbol,hs_meth_in_mouse_gbm$genes.Gene_Symbol)),y=genes_mouse_to_use),i)
  
list_datasets<-prepareDataTCGAandTurcan(hs_meth_in_mouse_gbm=hs_meth_in_mouse_gbm,hs_meth_in_mouse_lgg=hs_meth_in_mouse_lgg,genes_to_use=genes_mouse_to_use_rand,ann3=ann3,rand=T)

temp_tcga<-list_datasets[[1]]

tcga_dataset<-temp_tcga[which(temp_tcga$type%in%"LGG"),20:ncol(temp_tcga)]
tcga_dataset[,1]<-as.factor(tcga_dataset[,1])
colnames(tcga_dataset)[1]<-"status"
tcga_dataset[is.na(tcga_dataset)]<-0
tcga_dataset$status<-ifelse(tcga_dataset$status==1,"IDH_yes","IDH_no")

mat_turcan<-list_datasets[[2]]

common_genes_datasets<-intersect(colnames(tcga_dataset)[-1],colnames(mat_turcan)[-1])
tcga_dataset2<-tcga_dataset[,which(colnames(tcga_dataset)%in%c("status",common_genes_datasets))]
mat_turcan2<-mat_turcan[,which(colnames(mat_turcan)%in%c("status",common_genes_datasets))]

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/public_data")

#
# Use TCGA to build a model to predict Turcan
#

#
# Remove batch - effects
#

all_status<-rbind.fill(tcga_dataset2,mat_turcan2)[,1]
meth_all<-t(as.matrix(rbind.fill(tcga_dataset2,mat_turcan2)[,-1]))

groups<-as.factor(c(rep("TCGA",nrow(tcga_dataset2[,-1])),rep("Turcan",nrow(mat_turcan2[,-1]))))

combat_meth = ComBat(dat=meth_all, batch=groups)

combat_meth2<-data.frame(Groups=all_status,t(combat_meth))

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

overall_per<-confusionMatrix(rf_pred, as.factor(test_tcga_combat$IDH_status),positive="IDH_yes")$overall
byclass_per<-confusionMatrix(rf_pred, as.factor(test_tcga_combat$IDH_status),positive="IDH_yes")$byClass

all_stats<-rbind(all_stats,overall_per)

print(all_stats)

}

colnames(all_stats)<-names(overall_per)
  
all_stats2<-data.frame(size_rand_genes,all_stats)
write.table(all_stats2,file="TCGA_to_TURCAN_Combat_RandomGenes.txt",row.names=F,quote=F,sep="\t")