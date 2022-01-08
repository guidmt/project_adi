library(data.table)
library(minfi)
library("readxl")
library(factoextra)
library(ggpubr)

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

order_classes<-paste("s",1:81,sep="")
  
#
# Exploratory analysis
#
tab_input3<-tab_input3[-which(tab_input3$gene_symbol==""),]

gs<-strsplit(as.character(tab_input3[,1]),split="\\;")
gs2<-lapply(gs,FUN=unique)

gs3<-unlist(lapply(gs2,FUN=function(X){paste0(X,collapse=";")}))
tab_input3$gene_symbol<-gs3

check_illumina<-functionCheckGeneMet(gs2,c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2]))

# #how to check:
# xs<-as.numeric(check_illumina[[1]][37764])
# ys<- as.numeric(unlist(strsplit(check_illumina[[2]][37764],split=";")))
# 
# ytest<-c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2])
# gs2[xs]
# tab_input3[xs,1]
# ytest[ys]

# PCA hyper and hypo regulated genes
# xtest<-tab_input3[check_illumina[[1]],1]
# ytest<-c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2])
# common_check<-setdiff(xtest,ytest)


#
#hyper hypo
#
hyper_hypo_mat_meth<-tab_input3[check_illumina[[1]],-1]

melt_up_down<-melt(hyper_hypo_mat_meth)
melt_up_down_ann<-merge(melt_up_down,annotation_final,by.x="variable",by.y="sampleID")

pdf("turcan_hyper_hypo_meth_signals.pdf")

p2<-ggviolin(melt_up_down_ann, "epigenetic_status", "value", fill = "epigenetic_status",
             palette = c("red3", "green3"),
             add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p2+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

melt_meth_trim<-melt_up_down_ann[melt_up_down_ann$value>quantile(melt_up_down_ann$value)[2],]

p3<-ggviolin(melt_meth_trim, "epigenetic_status", "value", fill = "epigenetic_status",
             palette = c("red3", "green3"),
             add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p3+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

dev.off()


res.pca <- prcomp(t(hyper_hypo_mat_meth), scale = TRUE)

groups<-annotation_final[match(order_classes,annotation_final$sampleID),"epigenetic_status"]

groups<-as.factor(groups)

pdf("hyper_and_hypo_DMRs_Turcan.pdf")
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             legend.title = "Groups",
             label="none"
)
dev.off()

library(data.table)

hyper_hypo_mat_meth_genelevel<-setDT(tab_input3[check_illumina[[1]],])
hyper_hypo_mat_meth_genelevel2<- hyper_hypo_mat_meth_genelevel[, lapply(.SD, mean), by = gene_symbol]
hyper_hypo_mat_meth_genelevel3<-hyper_hypo_mat_meth_genelevel2[,-1]

melt_up_down<-melt(hyper_hypo_mat_meth_genelevel3)
melt_up_down_ann<-merge(melt_up_down,annotation_final,by.x="variable",by.y="sampleID")

pdf("turcan_hyper_hypo_meth_signals_genelevel.pdf",width=4,height=4)

p2<-ggviolin(melt_up_down_ann, "epigenetic_status", "value", fill = "epigenetic_status",
         palette = c("red3", "green3"),
         add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p2+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

melt_meth_trim<-melt_up_down_ann[melt_up_down_ann$value>quantile(melt_up_down_ann$value)[2],]

p3<-ggviolin(melt_meth_trim, "epigenetic_status", "value", fill = "epigenetic_status",
             palette = c("red3", "green3"),
             add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p3+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

dev.off()

res.pca <- prcomp(t(hyper_hypo_mat_meth_genelevel3), scale = TRUE)

pdf("hyper_and_hypo_DMRs_Turcan_genelevel.pdf")
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             legend.title = "Groups",
             label="none"
)
dev.off()

#
#hyper
#

# PCA all hyper in patients
check_illumina_up<-functionCheckGeneMet(gs2,c(hs_dmrs_up_from_mm[,2]))

hyper_mat_meth<-tab_input3[check_illumina_up[[1]],-1]

melt_up<-melt(hyper_mat_meth)
melt_up_ann<-merge(melt_up,annotation_final,by.x="variable",by.y="sampleID")

pdf("turcan_hyper_meth_signals.pdf")

p2<-ggviolin(melt_up_ann, "epigenetic_status", "value", fill = "epigenetic_status",
             palette = c("red3", "green3"),
             add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p2+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

melt_meth_trim<-melt_up_ann[melt_up_ann$value>quantile(melt_up_ann$value)[2],]

p3<-ggviolin(melt_meth_trim, "epigenetic_status", "value", fill = "epigenetic_status",
             palette = c("red3", "green3"),
             add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p3+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

dev.off()

res.pca <- prcomp(t(hyper_mat_meth), scale = TRUE)

groups<-annotation_final[match(order_classes,annotation_final$sampleID),"epigenetic_status"]

groups<-as.factor(groups)

pdf("hyper_DMRs_Turcan.pdf")
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             legend.title = "Groups",
             label="none"
)
dev.off()


hyper_mat_meth_genelevel<-setDT(tab_input3[check_illumina_up[[1]],])
hyper_mat_meth_genelevel2<- hyper_mat_meth_genelevel[, lapply(.SD, mean), by = gene_symbol]
hyper_mat_meth_genelevel3<-hyper_mat_meth_genelevel2[,-1]

melt_up<-melt(hyper_mat_meth_genelevel3)
melt_up_ann<-merge(melt_up,annotation_final,by.x="variable",by.y="sampleID")

pdf("turcan_hyper_meth_signals_genelevel.pdf")

p2<-ggviolin(melt_up_ann, "epigenetic_status", "value", fill = "epigenetic_status",
             palette = c("red3", "green3"),
             add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p2+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

melt_meth_trim<-melt_up_ann[melt_up_ann$value>quantile(melt_up_ann$value)[2],]

p3<-ggviolin(melt_meth_trim, "epigenetic_status", "value", fill = "epigenetic_status",
             palette = c("red3", "green3"),
             add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p3+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))


dev.off()

res.pca <- prcomp(t(hyper_mat_meth_genelevel3), scale = TRUE)

pdf("hyper_DMRs_Turcan_genelevel.pdf")
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             legend.title = "Groups",
             label="none"
)
dev.off()

#
#hypo
#

# PCA all hypo in patients
check_illumina_dw<-functionCheckGeneMet(gs2,c(hs_dmrs_down_from_mm[,2]))

hypo_mat_meth<-tab_input3[check_illumina_dw[[1]],-1]

melt_hypo<-melt(hypo_mat_meth)
melt_hypo_ann<-merge(melt_hypo,annotation_final,by.x="variable",by.y="sampleID")

pdf("turcan_hypo_meth_signals.pdf")

p2<-ggviolin(melt_hypo_ann, "epigenetic_status", "value", fill = "epigenetic_status",
             palette = c("red3", "green3"),
             add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p2+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

melt_meth_trim<-melt_hypo_ann[melt_hypo_ann$value>quantile(melt_hypo_ann$value)[2],]

p3<-ggviolin(melt_meth_trim, "epigenetic_status", "value", fill = "epigenetic_status",
             palette = c("red3", "green3"),
             add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p3+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

dev.off()

res.pca <- prcomp(t(hypo_mat_meth), scale = TRUE)

groups<-as.factor(groups)

pdf("hypo_DMRs_Turcan.pdf")
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             legend.title = "Groups",
             label="none"
)
dev.off()


hypo_mat_meth_genelevel<-setDT(tab_input3[check_illumina_dw[[1]],])
hypo_mat_meth_genelevel2<- hypo_mat_meth_genelevel[, lapply(.SD, mean), by = gene_symbol]
hypo_mat_meth_genelevel3<-hypo_mat_meth_genelevel2[,-1]

melt_hypo<-melt(hypo_mat_meth_genelevel3)
melt_hypo_ann<-merge(melt_hypo,annotation_final,by.x="variable",by.y="sampleID")

pdf("turcan_hypo_meth_signals_genesignal.pdf")

p2<-ggviolin(melt_hypo_ann, "epigenetic_status", "value", fill = "epigenetic_status",
             palette = c("red3", "green3"),
             add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p2+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

melt_meth_trim<-melt_hypo_ann[melt_hypo_ann$value>quantile(melt_hypo_ann$value)[2],]

p3<-ggviolin(melt_meth_trim, "epigenetic_status", "value", fill = "epigenetic_status",
             palette = c("red3", "green3"),
             add = "boxplot", add.params = list(fill = "white"))

my_comparisons<-list(c("CIMP+","CIMP-"))

print(p3+stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

dev.off()

res.pca <- prcomp(t(hypo_mat_meth_genelevel3), scale = TRUE)

pdf("hypo_DMRs_Turcan_genelevel.pdf")
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             legend.title = "Groups",
             label="none"
)
dev.off()

#
# Let's start with the machine learning algorithm
#

input_ml_adi_genes<-setDT(tab_input3[check_illumina[[1]],])
input_ml_adi_genes2<- input_ml_adi_genes[, lapply(.SD, mean), by = gene_symbol]
input_ml_adi_genes3<-input_ml_adi_genes2[,-1]

input_ml_adi_genes4<-t(input_ml_adi_genes3)
input_ml_adi_genes5<-data.frame(sampleID=rownames(input_ml_adi_genes4),input_ml_adi_genes4)
input_ml_adi_genes6<-merge(annotation_final,input_ml_adi_genes5,by.x="sampleID",by.y="sampleID")

input_ml_adi_genes7<-input_ml_adi_genes6[,c(3,14:ncol(input_ml_adi_genes6))]

rm("input_ml_adi_genes")
rm("input_ml_adi_genes2")
rm("input_ml_adi_genes3")

#
# ML1: Do  random forest
#

input_ml_adi_genes7$epigenetic_status<-ifelse(input_ml_adi_genes7$epigenetic_status=="CIMP+","methylator_positive","methylator_negative")

classes<-input_ml_adi_genes7$epigenetic_status
colnames(input_ml_adi_genes7)[1]<-"status"

input_ml_adi_genes7$status<-factor(input_ml_adi_genes7$status,levels=c("methylator_negative","methylator_positive"))

index <- createDataPartition(input_ml_adi_genes7[,"status"], p = 0.7, list = FALSE)

training<-input_ml_adi_genes7[index,]
test<-input_ml_adi_genes7[-index,]

library(ROCR)
library(caret)

print("trainControl")
ctrl <- trainControl(method="cv",
                     number=5,
                     classProbs=T,
                     savePredictions=T,
                     sampling = "down",
                     summaryFunction=twoClasSummary,
                     allowParallel=TRUE)
print("train")
rfmodel <- train(status ~ .,
                 data=training,
                 method="ranger",
                 trControl=ctrl,
                 num.trees=100)


rf_pred <- predict(rfmodel, test)

print("measure performances")

overall_per<-confusionMatrix(rf_pred, as.factor(test$status),positive="methylator_positive")$overall
byclass_per<-confusionMatrix(rf_pred, as.factor(test$status),positive="methylator_positive")$byClass

rf_pred <- predict(rfmodel, test,type="prob")

library(pROC)

rf.ROC <- roc(predictor=rf_pred$methylator_positive,response=test$status)
confidence_interval<-ci.auc(test$status, rf_pred$methylator_positive, conf.level = 0.9)
ci_string<-paste(round(as.numeric(confidence_interval),3)[2],"(",round(as.numeric(confidence_interval),3)[1],"-",round(as.numeric(confidence_interval),3)[3],")",sep="")

ciobj <- ci.se(rf.ROC, specificities=seq(0, 1, l=25), conf.level=0.95, boot.n=10000)

dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                     lower = ciobj[, 1],
                     upper = ciobj[, 3])

pdf("Turcan_DMRs_data_with_regions_Adi.pdf")

plot(rf.ROC,main=ci_string)
require(gridExtra)
per1<-tableGrob(data.frame(overall_per))
per2<-tableGrob(data.frame(byclass_per))
grid.arrange(per1,per2)

pAUCCI<-ggroc(rf.ROC) + theme_minimal() +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + coord_equal() +
  geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha= 0.2) + ggtitle(ci_string)

print(pAUCCI)

dev.off()
