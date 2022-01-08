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
# %of samples in which a gene is hyper or hypo methylated
# 
require("biomaRt")
#Warning-use the last version of Ensembl to convert mgi symbol to hg symbol - call useEnsembl to use a different version
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

dmrs_mm<-unique(dmrMouse[,19])
hs_dmrs_from_mm = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dmrs_mm , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

# get all the genes in the ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
all_hs_genes <- getBM(attributes=c("hgnc_symbol"),mart = mart)


# 
#  Compute % of samples with genes that are hyper and hypo-methylated in mouse and the two groups then you can create a boxplot of this
# 

#  Step 1: Download the methylation data in Human
current_dir_gbm<-paste("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/GDCdata","TCGA-GBM","harmonized/DNA_Methylation/Methylation_Beta_Value",sep="/")
mmgenes<-unique(hs_dmrs_from_mm[,2])
hs_meth_in_mouse_gbm<-create_matrix_meth(current_dir_gbm,genes_to_use=setdiff(all_hs_genes[,1],mmgenes))
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")
save(hs_meth_in_mouse_gbm,file="GBM_Mouse_beta_matrix_adi_pval0.05_lfc1.RANDOMGENES.RData")

current_dir_lgg<-paste("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/GDCdata","TCGA-LGG","harmonized/DNA_Methylation/Methylation_Beta_Value",sep="/")
hs_meth_in_mouse_lgg<-create_matrix_meth(current_dir_lgg,genes_to_use=setdiff(all_hs_genes[,1],mmgenes))
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")
save(hs_meth_in_mouse_lgg,file="LGG_Mouse_beta_matrix_adi_pval0.05_lfc1.RANDOMGENES.RData")
