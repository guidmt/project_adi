library(data.table)
library(stringr)
library(doParallel)
library(foreach)
require(gridExtra)
require(grid)
library(viridis)
library(ggridges)

source("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/create_matrix_meth2.R")
source("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/clean_methylation_data.R")

# 
#  Script that take the DMRs genes and compute the % of samples with hyper or hypo methylation across different TCGA-subgroups
# 

dmrMouse<-read.delim(file="/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/annotation_dmrs_adi_pval0.05_lfc1.txt",stringsAsFactors = F)

dmrMouse_rid<-unique(dmrMouse[,c(1,2,3,19,17)])
dmrMouse_stats<-dmrMouse[,c(19,13,14,17,18)]

# 
#  Extract the beta values in mouse 
# 

library(data.table)

setwd("/mnt/data/lab/gmt_data/data_brain/adi/CytosineReports")

list_meth_file<-grep(dir(),pattern="gz",value=T)

list_samples<-vector(mode='list',length(list_meth_file))

for(i in 1:length(list_meth_file)){
  
  print(i)
  
  print("reading data")
  meth_table_current_mouse<-fread(file=list_meth_file[i])
  colnames(meth_table_current_mouse)<-c("chromosome","position","strand","count_methylated","count_unmethylated","c_context","trinucleotide_context")
  
  beta_current_values<-data.frame()
  
  #select only the sites in which the gens were identified DEGs
  for(row_ann in 1:nrow(dmrMouse_rid)){
    
    print(paste("Genes:",row_ann))
    
    chrom<-dmrMouse_rid[row_ann,1]
    start<-dmrMouse_rid[row_ann,2]
    end<-dmrMouse_rid[row_ann,3]
    genes<-dmrMouse_rid[row_ann,4]
    
    temp_chr<-meth_table_current_mouse[chromosome==chrom]
    
    test<-data.frame(temp_chr[ position >=start & position <= end ])
    
    if(nrow(test)>=1){
      
      coordinates<-paste(range(test[,1]),collapse="-")
      # https://www.biostars.org/p/445286/
      beta_values<-sum(test$count_methylated)/(sum(test$count_methylated)+sum(test$count_unmethylated))
      
      df_beta<-data.frame(gene_symbol=genes,beta_values)
      
      beta_current_values<-rbind(beta_current_values,df_beta)
      
    }
    
  }
  
  
  list_samples[[i]]<-beta_current_values
  
}

res_meth<-do.call(cbind,list_samples[1:8])

res_meth2<-res_meth[,grep(colnames(res_meth),pattern="beta_values")]
colnames(res_meth2)<-sapply(strsplit(list_meth_file,split="\\."),"[[",1)

res_meth3<-as.data.frame(unique(cbind(genes=as.character(res_meth[,1]),res_meth2)))

remove_inf<-function(x){
  idx<-which(!is.finite(x))
  x[idx]<-0
  return(x)
}

res_meth3[,-1]<-sapply(res_meth3[,-1],remove_inf)

colnames(res_meth3)[1]<-"gene_symbol"

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")

save(res_meth3,file="Mouse_beta_matrix_adi_pval0.05_lfc1.RData")

# 
# %of samples in which a gene is hyper or hypo methylated
# 
require("biomaRt")
#Warning-use the last version of Ensembl to convert mgi symbol to hg symbol - call useEnsembl to use a different version
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

dmrs_up_mm<-unique(dmrMouse[dmrMouse$status=="DMRs-up",19])
hs_dmrs_up_from_mm = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dmrs_up_mm , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

dmrs_down_mm<-unique(dmrMouse[dmrMouse$status=="DMRs-down",19])
hs_dmrs_down_from_mm = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dmrs_down_mm , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

res_meth4<-res_meth3
res_meth4$status<-rep("NA",nrow(res_meth3))
res_meth4[res_meth4[,1]%in%dmrs_up_mm,"status"]<-"up"
res_meth4[res_meth4[,1]%in%dmrs_down_mm,"status"]<-"down"
# input_for_boxplot<-(melt(res_meth4))
# # ggboxplot(input_for_boxplot, "variable", "value", palette = c("#00AFBB", "#E7B800"),facet.by="status")
# 


# 
#  Compute % of samples with genes that are hyper and hypo-methylated in mouse and the two groups then you can create a boxplot of this
# 

#  Step 1: Download the methylation data in Human
current_dir_gbm<-paste("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/GDCdata","TCGA-GBM","harmonized/DNA_Methylation/Methylation_Beta_Value",sep="/")
hs_meth_in_mouse_gbm<-create_matrix_meth(current_dir_gbm,genes_to_use=c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2]))
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")
save(hs_meth_in_mouse_gbm,file="GBM_Mouse_beta_matrix_adi_pval0.05_lfc1.RData")

current_dir_lgg<-paste("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/GDCdata","TCGA-LGG","harmonized/DNA_Methylation/Methylation_Beta_Value",sep="/")
hs_meth_in_mouse_lgg<-create_matrix_meth(current_dir_lgg,genes_to_use=c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2]))
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")
save(hs_meth_in_mouse_lgg,file="LGG_Mouse_beta_matrix_adi_pval0.05_lfc1.RData")

#  Step 2: Filtering of methylation data in human
hs_meth_in_mouse_gbm_filt<-clean_methylation_data(mat_meth=hs_meth_in_mouse_gbm,genes_to_use_meth=c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2]))

hs_meth_in_mouse_lgg_filt<-clean_methylation_data(mat_meth=hs_meth_in_mouse_lgg,genes_to_use_meth=c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2]))

#  Step 3: Load annotation data methylation - define some groups
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human")

#load annotation data for TCGA
load("annotation_gbm_lgg.RData")

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis")
load("annotation_gbm_lgg_withmethylation.RData")

ann_meth<- ann_meth[,c("Samples","type","IDH1_R132H","TP53","ATRX","PDGFRA","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters","DNA_methylation_clusters","IDH_methylation_clusters","Supervised_DNA_methylation_clusters")]
ann_meth$IDH1_status<-rep("IDHR123H_negative",nrow(ann_meth))
ann_meth$IDH1_status[which(ann_meth$IDH1_R132H==1 & ann_meth$TP53==1 & ann_meth$ATRX==1)]<-"IDHR123H_positive"

ann_tot<-merge(ann_meth[,colnames(ann_meth)%in%c("Samples","IDH1_status","DNA_methylation_clusters","IDH_methylation_clusters","Supervised_DNA_methylation_clusters")],ann3,by.x="Samples",by.y="Samples")

#  Step 4: Combine all methylation
combined_meth<-combineMethMatrices(lgg_matrix=hs_meth_in_mouse_lgg_filt,gbm_matrix=hs_meth_in_mouse_gbm_filt)
colnames(hs_meth_in_mouse_lgg_filt)<-gsub(unlist(lapply(strsplit(colnames(hs_meth_in_mouse_lgg_filt),split="\\."),FUN=function(x){paste(x[1:3],collapse="-")})),pattern="\\.",replacement="-")
colnames(hs_meth_in_mouse_gbm_filt)<-gsub(unlist(lapply(strsplit(colnames(hs_meth_in_mouse_gbm_filt),split="\\."),FUN=function(x){paste(x[1:3],collapse="-")})),pattern="\\.",replacement="-")

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")
save(combined_meth,file="Combined_GBM_LGG_Mouse_beta_matrix_adi_pval0.05_lfc1.RData")

#
# Plot1: boxplot of the %of patients with a methylated regions
#

col_to_use <- c("IDH1_status","clusters","DNA_methylation_clusters","Supervised_DNA_methylation_clusters","IDH_methylation_clusters")
genes_to_consider_up <- hs_dmrs_up_from_mm[,2]
genes_to_consider_down <- hs_dmrs_down_from_mm[,2]

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")

plot_regions(combined_meth=combined_meth,ann_tot=ann_tot,genes_to_consider_up,genes_to_consider_down,col_to_use,thr=0.5,output_pdf="GBM_LGG_hyper_and_hypo_boxplot_on_tcga_filt_wilcox_0.05.pdf")
ann_lgg<-ann_tot[ann_tot$type%in%"LGG",]
plot_regions(combined_meth=hs_meth_in_mouse_lgg_filt[,c(4,8:ncol(hs_meth_in_mouse_lgg_filt))],ann_tot=ann_lgg,genes_to_consider_up,genes_to_consider_down,col_to_use,thr=0.5,output_pdf="LGG_hyper_and_hypo_boxplot_on_tcga_filt_wilcox_0.05.pdf")
ann_gbm<-ann_tot[ann_tot$type%in%"GBM",]
plot_regions(combined_meth=hs_meth_in_mouse_gbm_filt[,c(4,8:ncol(hs_meth_in_mouse_gbm_filt))],ann_tot=ann_gbm,genes_to_consider_up,genes_to_consider_down,col_to_use,thr=0.5,output_pdf="GBM_hyper_and_hypo_boxplot_on_tcga_filt_wilcox_0.05.pdf")

hs_meth_in_mouse_lgg_filt2<-hs_meth_in_mouse_lgg_filt[,c(4,8:ncol(hs_meth_in_mouse_lgg_filt))]
hs_meth_in_mouse_gbm_filt2<-hs_meth_in_mouse_gbm_filt[,c(4,8:ncol(hs_meth_in_mouse_gbm_filt))]

#
# Check if the DMRs in mouse are also DMRs in human
#


ann_lgg<-ann_tot[ann_tot$type%in%"LGG",]
ann_gbm<-ann_tot[ann_tot$type%in%"GBM",]

list_matrix_for_DMRs<-c("combined_meth","hs_meth_in_mouse_lgg_filt2","hs_meth_in_mouse_gbm_filt2")
list_ann<-c("ann_tot","ann_lgg","ann_gbm")
list_output<-c("GBM_LGG","LGG","GBM")
list_suffix<-c("R132H_TP53_ATRX","R132H_TP53_ATRX","R132H_TP53_ATRX")

for(my_analysis in 1:length(list_matrix_for_DMRs)){

    temp_ann<-get(list_ann[my_analysis])
    temp_matrix_analysis<-get(list_matrix_for_DMRs[my_analysis])
    current_prefix<-list_output[my_analysis]
    current_suffix<-list_suffix[my_analysis]

    samples_with_mutations<-temp_ann[temp_ann$IDH1_status=="IDHR123H_positive",1]
    samples_without_mutations<-temp_ann[temp_ann$IDH1_status=="IDHR123H_negative",1]
    
    if(length(samples_with_mutations) >10){
      
    df_stats2<-findHumanDMRs(temp_matrix_analysis,samples_with_mutations,samples_without_mutations)
    
    dmr_in_human<-df_stats2[which(df_stats2$padjust<=0.05),]
    dmr_in_human<-dmr_in_human[which(abs(dmr_in_human$all_delta)>=0.1),]
    
    dmr_in_human[,1]<-as.character( dmr_in_human[,1])
    
    library(VennDetail)
    
    up_in_human<-dmr_in_human[which(dmr_in_human$all_delta>= 0.1),1]
    down_in_human<-dmr_in_human[which(dmr_in_human$all_delta <= -0.1),1]
    
    ven <- venndetail(list(hyper_human = unique(up_in_human), 
                           hypo_human = unique(down_in_human),
                           hyper_mouse = unique(hs_dmrs_up_from_mm[,2]),
                           hypo_mouse = unique(hs_dmrs_down_from_mm[,2])))
    
    setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")
    output_file_pdf<-paste(paste(current_prefix,current_suffix,sep="_"),"Venny","pdf",sep=".")
    
    pdf(output_file_pdf)
    p1<-plot(ven)
    print(p1)
    p2<-plot(ven,type = "upset")
    print(p2)
    dev.off()
    
    wide_results<-result(ven, wide = TRUE)
    output_file_txt<-paste(paste(current_prefix,current_suffix,sep="_"),"Venny","txt",sep=".")
    write.table(wide_results,file=output_file_txt,sep="\t",row.names=F,quote=F)
    
    hyper_common<-wide_results[which(wide_results$hyper_human==1 & wide_results$hyper_mouse==1 & wide_results$SharedSets==2),1]     
    hypo_common<-wide_results[which(wide_results$hypo_human==1 & wide_results$hypo_mouse==1 & wide_results$SharedSets==2),1]     
    
    setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")
    col_to_use <- c("IDH1_status","clusters","DNA_methylation_clusters","Supervised_DNA_methylation_clusters","IDH_methylation_clusters")
    
    output_file2<-paste(paste(current_prefix,current_suffix,sep="_"),"TCGA_Clusters","pdf",sep=".")
    
    plot_regions(combined_meth=temp_matrix_analysis,ann_tot=temp_ann,hyper_common,hypo_common,col_to_use,thr=0.5,output_pdf=output_file2)
    
    }
    
}

#
# Do Principal Component Analysis of the mouse samples
#

melt_dmr<-melt(res_meth3[,-1])
colnames(melt_dmr)[2]<-"beta_values"

pdf("test.pdf")
ggplot(melt_dmr, aes(x = beta_values, y = variable, fill = stat(x))) +
  geom_density_ridges_gradient() +
  scale_fill_viridis_c(name = "Beta Values", option = "C")
dev.off()


library(factoextra)
library("FactoMineR")
library(ggfortify)
library(M3C)

input_pca2<-data.frame(groups=c("wt","wt","wt","mut","mut","mut"),t(unique(dmrMouse[,c(7:12)])))
colnames(input_pca2)[2:ncol(input_pca2)]<- paste("meth",seq(1,ncol(input_pca2)-1),sep="_")
  
input_pca2[,1]<-as.factor(input_pca2[,1])

respca2 <- prcomp(input_pca2[,-1], scale = TRUE)

pdf("mouse_prcomp_hyper_hypo.pdf")
autoplot(respca2, data = input_pca2, colour = 'groups', label = TRUE, label.size = 3)+ geom_hline(yintercept=0, linetype="dashed", color = "black")+ geom_vline(xintercept=0, linetype="dashed", color = "black")
dev.off()

pdf("PCA_GBM_LGG_hyper_hypo.pdf")
combined_meth_for_PCA<-t(combined_meth[,-1])
combined_meth_for_PCA2<-data.frame(ID=rownames(combined_meth_for_PCA),combined_meth_for_PCA)
combined_meth_for_PCA3<-merge(ann_tot,combined_meth_for_PCA2,by.x="Samples",by.y="ID",all.x=F,all.y=F)

temp_input_pca<-combined_meth_for_PCA3[,-c(1:20)]
temp_input_pca[is.na(temp_input_pca)]<-0

respca2 <- prcomp(temp_input_pca, scale = TRUE)

list_ann<-c("DNA_methylation_clusters","IDH_methylation_clusters","Supervised_DNA_methylation_clusters","IDH1_status","IDH1_R132H","type")
  
for(la in list_ann){
  
p<-autoplot(respca2, data = combined_meth_for_PCA3, colour = la, label = FALSE)+
         geom_hline(yintercept=0, linetype="dashed", color = "black")+ 
         geom_vline(xintercept=0, linetype="dashed", color = "black")
print(p+labs(title=la))

}

dev.off()


pdf("PCA_LGG_hyper_hypo.pdf")
combined_meth_for_PCA<-t(combined_meth[,-1])
combined_meth_for_PCA2<-data.frame(ID=rownames(combined_meth_for_PCA),combined_meth_for_PCA)
combined_meth_for_PCA3<-merge(ann_tot,combined_meth_for_PCA2,by.x="Samples",by.y="ID",all.x=F,all.y=F)
combined_meth_for_PCA3<-combined_meth_for_PCA3[combined_meth_for_PCA3$type=="LGG",]

temp_input_pca<-combined_meth_for_PCA3[,-c(1:20)]
temp_input_pca[is.na(temp_input_pca)]<-0

respca2 <- prcomp(temp_input_pca, scale = TRUE)

list_ann<-c("DNA_methylation_clusters","IDH_methylation_clusters","Supervised_DNA_methylation_clusters","IDH1_status","IDH1_R132H","type")

for(la in list_ann){
  
  p<-autoplot(respca2, data = combined_meth_for_PCA3, colour = la, label = FALSE)+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+ 
    geom_vline(xintercept=0, linetype="dashed", color = "black")
  print(p+labs(title=la))
  
}

dev.off()

#
# Approach similar to the paper between dogs and human
#
require(caret)
require(LibLinear)
require(LiblineaR)

#Measure the average Beta values for the genes in mus musculus
meth_data_mouse<-setDT(res_meth3)
meth_data_mouse2 <- data.frame(meth_data_mouse[, lapply(.SD, mean), by = gene_symbol])

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes_mm_map_hs = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = meth_data_mouse2[,1] , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

#Measure the average Beta values for the genes in homo sapiens
colnames(combined_meth)[1]<-"gene_symbol"
meth_data_hs<-setDT(combined_meth)
meth_data_hs2 <- data.frame(meth_data_hs[, lapply(.SD, mean), by = gene_symbol])

common_genes_hs_mm<-intersect(meth_data_hs2[,1],genes_mm_map_hs[,2])

mouse_genes_to_consider<-genes_mm_map_hs[genes_mm_map_hs[,2]%in%common_genes_hs_mm,1]
human_genes_to_consider<-genes_mm_map_hs[genes_mm_map_hs[,2]%in%common_genes_hs_mm,2]

meth_data_mouse3<-merge(genes_mm_map_hs,meth_data_mouse2[meth_data_mouse2[,1]%in%mouse_genes_to_consider,],by.x="MGI.symbol",by.y="gene_symbol")
meth_data_hs3<-meth_data_hs2[meth_data_hs2[,1]%in%human_genes_to_consider,]

# Now merge the two data.frame because is the only way to get the same number of columns
# Warning: i will not use more the mouse gene symbol but only the human gene symbol to be consistent

merged_matrix<-merge(meth_data_mouse3,meth_data_hs3,by.x="HGNC.symbol",by.y="gene_symbol")
# merged_matrix[,-c(1:2)]<-apply(merged_matrix[,-c(1:2)],1,scale)
# merged_matrix[is.na(merged_matrix)]<-0

idx1<-grep(colnames(merged_matrix),pattern="TCGA")
idx2<-grep(colnames(merged_matrix),pattern="HGNC.symbol")
human_tab<-merged_matrix[,c(idx2,idx1)]
colnames(human_tab)<-gsub(colnames(human_tab),pattern="\\.",replacement="-")
                          
input_for_ml<-t(human_tab[,-1])
colnames(input_for_ml)<-human_tab[,1]
input_for_ml2<-data.frame(ID=rownames(input_for_ml),input_for_ml)
input_for_ml3<-merge(ann_tot,input_for_ml2,by.x="Samples",by.y="ID",all.x=F,all.y=F)
colnames(input_for_ml3)<-gsub(colnames(input_for_ml3),pattern="\\.",replacement="-")

idx1<-grep(colnames(merged_matrix),pattern="TCGA",invert=T)
idx2<-grep(colnames(merged_matrix),pattern="HGNC.symbol")
mouse_tab<-merged_matrix[,c(idx2,idx1)][,-c(2:4)] #remove columns with gene-symbol and control EC

input_mm<-t(mouse_tab[,-1])
colnames(input_mm)<-mouse_tab[,1]
input_mm2<-data.frame(ID=c("mut","mut","mut","wt","wt","wt","wt"),input_mm)
colnames(input_mm2)[2:ncol(input_mm2)]<-gsub(colnames(input_mm2)[2:ncol(input_mm2)],pattern="\\.",replacement="-")
colnames(input_mm2)[1]<-"status_mouse"
rownames(input_mm2)<-colnames(mouse_tab)[-1]

list_ann<-c("DNA_methylation_clusters","IDH_methylation_clusters","Supervised_DNA_methylation_clusters","histology","IDH1_status","IDH1_R132H")

lgg_ml_res<-libLinear_mouse_hs(input_for_ml3,input_mm2,list_ann=list_ann,type_tumours="LGG",prefix_output="LGG",suffix_output="DMRs_hyper_hypo")
gbm_lgg_ml_res<-libLinear_mouse_hs(input_for_ml3,input_mm2,list_ann=list_ann,type_tumours=c("GBM","LGG"),prefix_output="GBM_LGG",suffix_output="DMRs_hyper_hypo")

lgg_temp<-lgg_ml_res[[2]]
names(lgg_temp)<-paste("ID",1:length(lgg_temp),sep="")
lgg_input_for_boxplot<-melt(lgg_temp)

pdf("LGG_boxplot_probabilities.pdf")
p1<-ggboxplot(lgg_input_for_boxplot, 
          "variable", "value", 
          orientation = "horizontal",
          ylab="Probabilities")+facet_grid(rows="L1",space="free",scales="free")
print(p1)
p2<-ggboxplot(lgg_input_for_boxplot, 
          "variable", "value", 
          orientation = "horizontal",
          ylab="Probabilities")
print(p2)

p3<-ggplot(lgg_input_for_boxplot, aes(x = value, y = variable, fill = stat(x))) +
  geom_density_ridges_gradient() +
  scale_fill_viridis_c(name = "Probabilities", option = "C")+geom_vline(xintercept = 0.5, linetype="dotted", 
                                                                        color = "blue", size=1)

print(p3)

dev.off()


gbm_lgg_temp<-gbm_lgg_ml_res[[2]]
names(gbm_lgg_temp)<-paste("ID",1:length(gbm_lgg_temp),sep="")
gbm_lgg_input_for_boxplot<-melt(gbm_lgg_temp)

# 
# input_for_density<-gbm_lgg_input_for_boxplot[gbm_lgg_input_for_boxplot[,4]%in%c("ID1","ID5","ID6"),]
# 

pdf("GBM_LGG_boxplot_probabilities.pdf")
p1<-ggboxplot(gbm_lgg_input_for_boxplot, 
              "variable", "value", 
              orientation = "horizontal",
              ylab="Probabilities")+facet_grid(rows="L1",space="free",scales="free")
print(p1)
p2<-ggboxplot(gbm_lgg_input_for_boxplot, 
              "variable", "value", 
              orientation = "horizontal",
              ylab="Probabilities")
print(p2)

# 
# input_for_density<-gbm_lgg_input_for_boxplot[gbm_lgg_input_for_boxplot[,4]%in%c("ID1","ID5","ID6"),]
# 
p3<-ggplot(gbm_lgg_input_for_boxplot, aes(x = value, y = variable, fill = stat(x))) +
  geom_density_ridges_gradient() +
  scale_fill_viridis_c(name = "Probabilities", option = "C")+geom_vline(xintercept = 0.5, linetype="dotted", 
                                                                        color = "blue", size=1)

print(p3)

dev.off()


save(list=c("lgg_ml_res","gbm_ml_res"),file="res_liblinear_ML_mouse_human.RData")

