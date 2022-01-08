library(data.table)
library(minfi)
library("readxl")
library(factoextra)

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

matrix_cell_line<-tab_input[,colnames(tab_input)%in%c("Gene Symbol",grep(colnames(tab_input),pattern="^P",value=T))]
matrix_cell_line2<-cbind(gene_symbol=matrix_cell_line[,ncol(matrix_cell_line)],matrix_cell_line[,-ncol(matrix_cell_line)])
matrix_cell_line2[,-1]<-apply(matrix_cell_line2[,-1],2,as.numeric)

#
# Do a PCA analysis of the hyper and hypo genes using the 81 samples
#

my_data <- data.frame(read_excel("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/public_data/Table 5a.xls"))
colnames(my_data)<-my_data[1,]
my_data2<-my_data[-1,]
my_data2[,1]<-paste("s",my_data2[,1],sep="")
my_data2<-my_data2[1:81,]

# PCA all hyper in patients
hyper_mat_meth<-tab_input3[tab_input3[,1]%in%hs_dmrs_up_from_mm[,2],-1]

res.pca <- prcomp(t(hyper_mat_meth), scale = TRUE)

groups<-as.factor(my_data2[,2])

pdf("hyper_DMRs_Turcan.pdf")
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             legend.title = "Groups",
             label="none"
)
dev.off()

var_hyper<-apply(hyper_mat_meth,1,var)

what_select<-names(var_hyper[1:10000])

hyper_mat_meth_var<-hyper_mat_meth[rownames(hyper_mat_meth)%in%what_select,]

res.pca <- prcomp(t(hyper_mat_meth_var), scale = TRUE)

groups<-as.factor(my_data2[,2])

pdf("hyper_DMRs_Turcan_var.pdf")
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             legend.title = "Groups",
             label="none"
)
dev.off()


# PCA all hypo in patients
hypo_mat_meth<-tab_input3[tab_input3[,1]%in%hs_dmrs_down_from_mm[,2],-1]

res.pca <- prcomp(t(hypo_mat_meth), scale = TRUE)

groups<-as.factor(my_data2[,2])

pdf("hypo_DMRs_Turcan.pdf")
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             legend.title = "Groups",
             label="none"
)
dev.off()


var_hypo<-apply(hypo_mat_meth,1,var)

what_select<-names(var_hypo[1:1000])

hypo_mat_meth_var<-hypo_mat_meth[rownames(hypo_mat_meth)%in%what_select,]

res.pca <- prcomp(t(hypo_mat_meth_var), scale = TRUE)

groups<-as.factor(my_data2[,2])

pdf("hypo_DMRs_Turcan_var.pdf")
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             legend.title = "Groups",
             label="none"
)
dev.off()

