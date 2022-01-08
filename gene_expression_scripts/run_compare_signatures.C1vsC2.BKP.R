library(biomaRt)
library(sva)
library(edgeR)
library(pheatmap)
library(sva)
library(glmnet)
library(ComplexHeatmap)
library(data.table)
library(circlize)
library(stringr)
library(stringi)


system("cp run_compare_signatures.C1vsC2.R run_compare_signatures.C1vsC2.BKP.R")

source("convert_mouse_to_human.R")

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis")
ortholog_table<-read.delim(file='HMD_HumanPhenotype.rpt',header=F) #see https://www.sciencedirect.com/science/article/pii/S1074761319301268#fig2

setwd("/mnt/data/lab/gmt_data/data_brain/")
mm_sampleinfo<-read.delim(file='sampleinfo2.txt')

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis")
# Upload the RNA-seq values from Human
load('tcga.lgg.mut.RData')
load('tcga.gbm.mut.RData')

input_ge<-rbind(tcga.gbm,tcga.lgg)

col_to_select<-c(colnames(input_ge[,c(1:6)]),'Type')

setwd('/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human')

#load annotation data
load("annotation_gbm_lgg.RData")
#load annotation 
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/public_data")
load("CGGA_ann.RData")
#load microarray for metanalysis
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/public_data")
load("database_microarray_brain.RData")


basic_directory<-c("DESeq2_rmv5_CSHC1A_Run1C2_Adi18","DESeq2_rmv5_rmvAdi18","DESeq2_rmv5_rmvCSHC1A_rmvAdi18")
list_directories<-c("Comp1_C1_vs_C2","Comp1_C1_vs_C2","Comp1_C1_vs_C2")
comparisons<-c('C1:C2','C1:C2','C1:C2')

list_pval<-c(0.001,0.001,0.001,0.001)
list_lfc<-c(2,2,2,2)

for(i in 1:length(list_directories)){

	bd<-basic_directory[i]
	pval<-list_pval[i]
	lfc<-list_lfc[i]

	print(pval)
	#
	# Define where save the output file
	#

	setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_17_08_2020/res_17_08_2020")
	system(paste("mkdir",bd))

	setwd(paste("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_17_08_2020/res_17_08_2020",bd,sep="/"))
        output_dir<-getwd()
	#end

	#define the working directory
	current_dir<-paste("/mnt/data/lab/gmt_data/data_brain/ForAdi_20200727",bd,list_directories[i],sep='/')

	setwd(current_dir)

	input_file<-grep(dir(),pattern='.sigDE.csv',value=T)

	degs_in_mouse<-read.csv(file=input_file,sep=',')

	list_genes<-as.character(degs_in_mouse$Gene_name)

	#1) Get the groups in this comparison
	
	a<-unlist(strsplit(comparisons[i],split='\\:'))[1]
	b<-unlist(strsplit(comparisons[i],split='\\:'))[2]

	#2) Get the samples in this comparison

	statistics_degs<-degs_in_mouse[,c(1:8)]

	matrix_mm<-degs_in_mouse[,-c(1:8)]
	
	rownames(matrix_mm)<-paste(degs_in_mouse$ID,degs_in_mouse$Gene_name,sep=';')

	current_samples<- colnames(matrix_mm)
	
	mm_sampleinfo_rid<-mm_sampleinfo[which(mm_sampleinfo$Sample_name%in%current_samples),]
	conds<-unique(mm_sampleinfo_rid$Condition2)

	matrix_mm2<-data.frame(matrix_mm[,match(mm_sampleinfo_rid$Sample_name,colnames(matrix_mm))])

	ann_mm_clean<-mm_sampleinfo_rid[,c('Condition','Condition2')]

	rownames(ann_mm_clean)<-mm_sampleinfo_rid$Sample_name

	set.seed(123) 
	
	#	
	# 4) Filter based on DEGs analysis 
	#

	statistics_degs<-degs_in_mouse[,c(1:8)]
	genes_to_save<-grep(grep(grep(grep(grep(grep(statistics_degs$Gene_name,pattern='Rik',invert=T,value=T),pattern='Gm',value=T,invert=T),pattern='AC',invert=T,value=T),pattern='CT',invert=T,value=T),pattern="mt-",invert=T,value=T),pattern="CAA",invert=T,value=T)

	statistics_degs<-statistics_degs[which(statistics_degs$Gene_name%in%genes_to_save),]

	selected_genes_with_degs<-as.character(statistics_degs[statistics_degs$padj<=pval & abs(statistics_degs$log2FoldChange)>=lfc,1])

	up_genes<-statistics_degs[statistics_degs$padj<=pval & statistics_degs$log2FoldChange>=lfc,1]
	down_genes<-statistics_degs[statistics_degs$padj<=pval & statistics_degs$log2FoldChange <= -c(lfc),1]

	output_degs<-paste(paste(list_directories[i],pval,lfc,sep='_'),'.txt',sep='')

	output_txt<-paste(paste(output_dir,output_degs,sep="/"))
	
	write.table(selected_genes_with_degs,file=output_txt,sep='\t',row.names=T,quote=F)

	ortholog_table_mgi<-ortholog_table[which(ortholog_table[,5]%in%selected_genes_with_degs),]

	ortholog_table_ensembl<-unique(convertMouseGeneList(selected_genes_with_degs))

	ortholog_table_res<-integrated_conversion_mouse_human(mgi=ortholog_table_mgi,ensembl=ortholog_table_ensembl)

	output_degs2<-paste(paste(list_directories[i],pval,lfc,sep='_'),'.ort.txt',sep='')
	output_txt2<-paste(paste(output_dir,output_degs2,sep="/"))

	hs_genes_sig<-unique(ortholog_table_res[,1])
	mm_genes_sig<-unique(ortholog_table_res[,2])

	hs_genes_sig_up<-ortholog_table_res[ortholog_table_res[,2]%in%up_genes,1]
 	hs_genes_sig_down<-ortholog_table_res[ortholog_table_res[,2]%in%down_genes,1]

	ortholog_table_res$status[ortholog_table_res[,1]%in%hs_genes_sig_up]<-"up"
	ortholog_table_res$status[ortholog_table_res[,1]%in%hs_genes_sig_down]<-"down"
  	
	write.table(ortholog_table_res,file=output_txt2,sep='\t',row.names=T,quote=F)

  	#
  	# Create matrix h.sapiens --NOT-- processed with combat
	#
	col_to_select<-c(colnames(input_ge[,c(1:6)]),'Type')
  	matrix_to_plot_hs_nocombat<-t(input_ge[,-c(which(colnames(input_ge)%in%col_to_select))])
	colnames(matrix_to_plot_hs_nocombat)<-unlist(lapply(X=strsplit(input_ge$Sample,split='\\-'),FUN=function(X)paste(X[1:3],collapse='-')))

	matrix_to_plot_hs_nocombat2<-matrix_to_plot_hs_nocombat[which(rownames(matrix_to_plot_hs_nocombat)%in%hs_genes_sig),]
	
	#
	# Create matrix mouse with all genes derived from the DEGs analysis
	#

	matrix_mm3<-data.frame(gs=unlist(lapply(strsplit(rownames(matrix_mm2),split=';'),'[[',2)),matrix_mm2)

	matrix_to_plot_ms<-setDT(matrix_mm3[matrix_mm3$gs%in%selected_genes_with_degs,])

    	matrix_to_plot_ms <- as.data.frame(matrix_to_plot_ms[, lapply(.SD, mean), by = gs])

	rownames(matrix_to_plot_ms)<-matrix_to_plot_ms[,1]

	matrix_to_plot_ms<-matrix_to_plot_ms[,-1]
	#
	# Create a matrix with the genes that are in common with humans
	#

	matrix_to_plot_ms_ort<-setDT(matrix_mm3[matrix_mm3$gs%in%mm_genes_sig,])
        matrix_to_plot_ms_ort <- as.data.frame(matrix_to_plot_ms_ort[, lapply(.SD, mean), by = gs])

	rownames(matrix_to_plot_ms_ort)<-matrix_to_plot_ms_ort[,1]

	matrix_to_plot_ms_ort<-matrix_to_plot_ms_ort[,-1]
	
	setwd('/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_17_08_2020/')
	source("support_functions_heatmaps.R")
        source("support_functions_MutBurden_GE_METH.R")
        source("create_matrix_meth.R")
        source("support_function_heatmap_mouse.R")

	###############################################
	# Plot human results
	################################################ 
	
	setwd(output_dir)

	annotation_all<-ann3[,c("Samples","type","IDH1_R132H","TP53","ATRX","PDGFRA","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]

	# No Combat: all cancers
	output_hs<-paste(paste(list_directories[i],'method.fc',pval,lfc,sep='.'),'.nocombat.png',sep='')
  	output_hs_var<-paste(paste(list_directories[i],'method.fc',pval,lfc,sep='.'),'.nocombat.var.png',sep='')

	plot_heatmap_hs(matrix_to_plot_hs_nocombat2,annotation_all,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=TRUE,number_clusters_row=2)	

	library(ConsensusClusterPlus)
	results = ConsensusClusterPlus(matrix_to_plot_hs_nocombat2,
				       maxK=10,
				       reps=50,
				       pItem=0.8,
				       pFeature=1,
				       title="test",
				       clusterAlg="km",
				       seed=1262118388.71279,plot="png")

	classesdf<-data.frame(ID=names(results[[8]]$consensusClass),
			      classes=results[[8]]$consensusClass)
	
	library(glmnet)
	cvfit = cv.glmnet(t(matrix_to_plot_hs_nocombat2),classesdf$classes,type.measure="class",alpha=1,family="multinomial")
	test<-coef(cvfit, s = "lambda.min")

	resLasso<-function(x){
	idx<-which(x!=0)
	cg<-names(x[idx,])[-1]
	return(cg)
	}

	signatureGenes<-unique(unlist(lapply(test,resLasso)))

	output_hs<-paste(paste(list_directories[i],'method.fc',pval,lfc,sep='.'),'.nocombat.test.png',sep='')	
	output_hs_var<-paste(paste(list_directories[i],'method.fc',pval,lfc,sep='.'),'.nocombat.var.test.png',sep='')
		       
       	plot_heatmap_hs(matrix_to_plot_hs_nocombat2[rownames(matrix_to_plot_hs_nocombat2)%in%signatureGenes,],annotation_all,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,km=F,split=T,how_split=classesdf$classes,number_clusters_row=2)

	#
	# See the methylation data
	#

	run_meth=FALSE

	if(run_meth==TRUE){

	setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis")
	
	load("annotation_gbm_lgg_withmethylation.RData")
	        
	output_hs<-paste(paste(list_directories[i],'method.fc',pval,lfc,sep='.'),'.ge.meth.nocombat.png',sep='')

	output_hs_meth<-paste(paste(list_directories[i],'method.fc',pval,lfc,sep='.'),'.ONLY.METH.nocombat.png',sep='')

	ann_meth<- ann_meth[,c("Samples","type","IDH1_R132H","TP53","ATRX","PDGFRA","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters","IDH_methylation_clusters","Supervised_DNA_methylation_clusters")]

	if(length(grep(dir(),pattern=paste(list_directories[i],pval,lfc,"meth.RData",sep=".")))==0){

        # get the methylation values 
	print("create methylation matrix GBM")

	genes_to_use<-rownames(matrix_to_plot_hs_nocombat2)

	meth_gbm<-create_matrix_meth("TCGA-GBM",genes_to_use)

	print("create methylation matrix LGG")

	meth_lgg<-create_matrix_meth("TCGA-LGG",genes_to_use)

	mat_meth<-cbind(meth_gbm,meth_lgg[,-c(1:2)])
	colnames(mat_meth)<-gsub(unlist(lapply(strsplit(colnames(mat_meth),split="\\."),FUN=function(x){paste(x[1:3],collapse="-")})),pattern="\\.",replacement="-")
	save(mat_meth,file=paste(list_directories[i],pval,lfc,"meth.RData",sep="."))
	
	} else {
	
	load(paste(list_directories[i],pval,lfc,"meth.RData",sep="."))
	
	}
	setwd(output_dir)

	plot_heatmap_hs_ge_meth(matrix_to_plot_hs_nocombat2,
				mat_meth=mat_meth,
				ann_meth,
				output_hs=output_hs,
				output_hs2=output_hs_meth,
				column_title="DEGs genes and associated meth islands",
				width=15,
				height=10,
				split=T,
				km=TRUE,
				qthr=30,
				number_clusters=3,
				number_clusters_row=1)

	}

	##########################################

	
	###############################################
	# Plot mouse results: all degs genes
        ###############################################	

	setwd(output_dir)

	output_mm<-paste('Mus_musculus',list_directories[i],'.method.fc',pval,lfc,'.png',sep='')

	cal_z_score <- function(x){
		  (x - mean(x)) / sd(x)
	}

	mat_mm<-log(matrix_to_plot_ms+1,2)
	mat_mm<-t(apply(mat_mm,1,cal_z_score))
	colnames(mat_mm)<-colnames(matrix_to_plot_ms)
	
	colours<-list("Condition"=c("C1A"="orange","C1B"="blue","C2"="yellow"),
		   "Condition2"=c("C1"="red","C2"="green"))
	
	colAnn <- HeatmapAnnotation(df=ann_mm_clean, which="col",col=colours,
                                    annotation_width=unit(c(1, 4), "cm"),
                                    gap=unit(1, "mm"),
                                    annotation_name_side = "left",
                                    border=TRUE)


	png(output_mm, height=8,width=8,units='in',res=300)

	mmmap<-Heatmap(mat_mm,
		show_column_names=FALSE,
		show_row_names=F,
		cluster_columns = T,
		show_column_dend = TRUE,
                row_dend_reorder = F,
                column_dend_reorder = F,
                show_row_dend = TRUE,
                clustering_distance_rows = "euclidean",
		clustering_method_rows = "complete",
		clustering_distance_columns = "euclidean",
	        clustering_method_columns = "complete",
	        top_annotation=colAnn)

	draw(mmmap, heatmap_legend_side="bottom", annotation_legend_side="right")

	dev.off()
        
	###############################################
        # Plot mouse results: ort genes with human
        ############################################### 

	output_mm<-paste('Mus_musculus',list_directories[i],'.method.ort.fc',pval,lfc,'.png',sep='')

	cal_z_score <- function(x){
		  (x - mean(x)) / sd(x)
	}

	mat_mm<-log(matrix_to_plot_ms_ort+1,2)
	mat_mm<-t(apply(mat_mm,1,cal_z_score))
	colnames(mat_mm)<-colnames(matrix_to_plot_ms_ort)

	colours<-list("Condition"=c("C1A"="orange","C1B"="blue","C2"="yellow"),
		                    "Condition2"=c("C1"="red","C2"="green"))

	colAnn <- HeatmapAnnotation(df=ann_mm_clean, which="col",col=colours,
	                            annotation_width=unit(c(1, 4), "cm"),
	                            gap=unit(1, "mm"),
                                    annotation_name_side = "left",
			            border=TRUE)

	png(output_mm, height=8,width=8,units='in',res=300)

	mmmap<-Heatmap(mat_mm,
                       show_column_names=FALSE,
		       show_row_names=F,
		       cluster_columns = T,
		       show_column_dend = TRUE,
	               row_dend_reorder = F,
		       column_dend_reorder = F,
		       show_row_dend = TRUE,
		       clustering_distance_rows = "euclidean",
		       clustering_method_rows = "complete",
		       clustering_distance_columns = "euclidean",
		       clustering_method_columns = "complete",
		       top_annotation=colAnn)
        
	draw(mmmap, heatmap_legend_side="bottom", annotation_legend_side="right")

	dev.off()

	setwd("/mnt/data/lab/gmt_data/data_brain/ForAdi_20200727")
	  
}	
