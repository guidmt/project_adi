# load libraries to use during the analysis

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
library(xlsx)
library(TCGAbiolinks)
library(caret)
library(ROCR)
library(ranger)
library(pROC)
library(GSVA)
library(survival)
library(RegParallel)

# load functions useful for the analysis

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/")
source("convert_mouse_to_human.R")
source("support_functions_heatmaps.R")
source("support_functions_MutBurden_GE_METH.R")
source("create_matrix_meth.R")
source("support_function_heatmap_mouse.R")
source("heatmapMouse.R")
source("support_function_survival.R")
source("compareDegsHs_MM.R")
source("mlForModel.R")
source("CrossSpeciesComparison.R")

# Define a function for gsva analysis
funGSVA<-function(file_gsva,X,min.sz=1,max.sz=2000,cut_off_gsva=0.5,cut_off_logfc=1,method="singscore",cox_genes=F,list_genes_cox=NULL){
  
  genes_for_gsva<-read.delim(file=file_gsva,stringsAsFactors=F)
  
  ovp_gsva_fs<-length(intersect(genes_for_gsva[,1],rownames(X)))
  
  print(paste("The overlap between C1AC1B and X is:",ovp_gsva_fs))
  
  library(GSVA)
  library(singscore)
  
  if(cox_genes==F){
    
  up_genes<-genes_for_gsva[genes_for_gsva[,3]%in%"up",1]
  down_genes<-genes_for_gsva[genes_for_gsva[,3]%in%"down",1]
  
  }else{
    
  genes_for_gsva<-genes_for_gsva[genes_for_gsva[,1]%in%list_genes_cox,]
  up_genes<-genes_for_gsva[genes_for_gsva[,3]%in%"up",1]
  down_genes<-genes_for_gsva[genes_for_gsva[,3]%in%"down",1]
  
  }
  
  if(method=="gsva"){
    
    if(length(down_genes)!=0 & length(up_genes)!=0){
      
    gsva_up <- gsva(as.matrix(X), list(up_genes),min.sz=min.sz, max.sz=max.sz, verbose=TRUE)
    gsva_down<- gsva(as.matrix(X), list(down_genes),min.sz=min.sz, max.sz=max.sz, verbose=TRUE)
    
    scores_combined<-gsva_up-gsva_down
    
    cut_off <- cut_off_gsva
    
    } 
    
    if(length(down_genes)==0){
      
      # if there are only up regulated genes
      gsva_up <- gsva(as.matrix(X), list(up_genes),min.sz=min.sz, max.sz=max.sz, verbose=TRUE)
      # the scores combined is given by the scores of the up-regulated genes
      scores_combined<-gsva_up
      
      # the scores of the down regulated genes is equal to the scores of the up-regulated genes
      gsva_down<-t(as.matrix(rep(0,length(gsva_up))))

      cut_off <- cut_off_gsva
      
    }
    
    if(length(up_genes)==0){
      
      # if there are only up regulated genes
      gsva_down <- gsva(as.matrix(X), list(down_genes),min.sz=min.sz, max.sz=max.sz, verbose=TRUE)
      # the scores combined is given by the scores of the up-regulated genes
      scores_combined<-gsva_down
      
      # the scores of the down regulated genes is equal to the scores of the up-regulated genes
      gsva_up<-t(as.matrix(rep(0,length(gsva_down))))
      
      cut_off <- cut_off_gsva
      
    }
    
    status_C1AC1B<-rep("NA",length(scores_combined))
    
    status_C1AC1B[which(scores_combined > cut_off_gsva)]<-"high"
    status_C1AC1B[which(scores_combined <= cut_off_gsva)]<-"low"    

    names(status_C1AC1B)<-colnames(X)
    
    output_df<- data.frame(ID=colnames(X),status_C1AC1B,up_score=t(gsva_up),down_score=t(gsva_down),global_score=t(scores_combined))
    
  
    colnames(output_df)[2]<-"status_C1AC1B"

  }

  if(method=="singscore"){  
  
    xdt<-setDT(data.frame(Gene=rownames(X),X))
    xdt_mean <- xdt[, lapply(.SD, mean), by = Gene]
    xdt_mean_df <- data.frame(xdt_mean)
    rownames(xdt_mean_df)<-xdt_mean_df[,1]
    rankData <- rankGenes(xdt_mean_df[,-1])
    scoredf <- simpleScore(rankData, upSet = unique(up_genes), downSet = unique(down_genes))
    
    p<-plotRankDensity(rankData[,2,drop = FALSE], upSet = unique(up_genes), downSet = unique(down_genes), isInteractive = FALSE)
    
    print(p)
    
    output_df<-data.frame(ID=colnames(X),
                          status_C1AC1B="",
                          up_score=scoredf$UpScore,
                          down_score=scoredf$DownScore,
                          global_score=scoredf$TotalScore)
    
    colnames(output_df)[2]<-"status_C1AC1B"

  }

  if(method=="logFC"){
  
  mean_up<-apply(X[rownames(X)%in%up_genes,],2,mean)
  
  #it there are genes down-regulated
  if(length(down_genes)!=0){
    
  # and if there are more than 2 down-regulated genes use apply
      if(length(down_genes)>=2){
      
      mean_down<-apply(X[rownames(X)%in%down_genes,],2,mean)
    
      } else {
  # otherwise use the mean function
      mean_down<-mean(X[rownames(X)%in%down_genes,])
      
      }
  
  } else {
    
    # it there are not down regulated genes the score is 0
    mean_down<-t(as.matrix(rep(0,length(mean_up))))

  }
  
  if(length(down_genes)!=0){
    
  logFC<-mean_up-mean_down
  
  } else {
    
  logFC<-mean_up
  
  }
  
  status_C1AC1B<-rep("NA",length(logFC))
  
  cut_off <- cut_off_logfc
  
  status_C1AC1B[which(logFC >= cut_off)]<-"high"    
  status_C1AC1B[which(logFC < cut_off)]<-"low"

  output_df<- data.frame(ID=colnames(X),
			 status_C1AC1B=status_C1AC1B,
			 up_score=as.numeric(mean_up),
			 down_score=as.numeric(mean_down),
			 global_score=as.numeric(logFC))
  
  colnames(output_df)[2]<-"status_C1AC1B"

  }

  output_list<-vector(mode="list",2)
  output_list[[1]]<-paste(method,"MCP",cut_off,sep="_")
  output_list[[2]]<-output_df

  return(output_list)
  
}



#
# Run a code to debug the data
#

x<-date()
xdate<-paste("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/bkp_script/","run_compare_signatures.C1vsC2.GENERAL.",unlist(lapply(strsplit(paste(x),split=" "),FUN=function(x){paste(x[1:6],collapse="_")})),".R",sep="")
system(paste("cp","/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/run_compare_signatures.GENERAL.R",xdate))

#bkp all script
system(paste("cp","*.R","/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/bkp_script/"))


# Upload a table with the orthologs genes 
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis")
ortholog_table<-read.delim(file="HMD_HumanPhenotype.rpt",header=F) #see https://www.sciencedirect.com/science/article/pii/S1074761319301268#fig2

# Upload a sampleinfo with the mouse samples
mm_sampleinfo<-read.delim(file="/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/sampleinfo_04_02_2021.txt",check.names=FALSE)

# Upload the RNA-seq values from Human
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis")
load("tcga.lgg.mut.RData")
load("tcga.gbm.mut.RData")

input_ge<-rbind(tcga.gbm,tcga.lgg)
col_to_select<-c(colnames(input_ge[,c(1:6)]),"Type")

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human")

#load annotation data for TCGA
load("annotation_gbm_lgg_July2021.RData")
colnames(ann3)[1]<-"Samples"

# read the sampleinfo to define which analysis perform 
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021")

sampleinfo<-read.delim(file="sampleinfo_C1vsC2.txt",header=F,stringsAsFactors=F,comment.char="#")

for(i in 1:nrow(sampleinfo)){
	
	bd<-sampleinfo[i,1]
	directory_exp<-sampleinfo[i,2]
	pval<-sampleinfo[i,4]
	lfc<-sampleinfo[i,5]
	mkcluster<-sampleinfo[i,6]
	feature_selection<-sampleinfo[i,7]
	do_gsva<-sampleinfo[i,8]

	print(pval)
	
	#
	# Define where save the output file
	#

	setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/res_04_02_2021")
	system(paste("mkdir -p",paste(bd,directory_exp,sep="/")))

	setwd(paste("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/res_04_02_2021",paste(bd,directory_exp,sep="/"),sep="/"))
	output_dir<-getwd()
	
	#
	#define the working directory
	#
	
	current_dir<-paste("/mnt/data/lab/gmt_data/data_brain",bd,directory_exp,sep="/")

	setwd(current_dir)

	input_file<-grep(dir(),pattern=".sigDE.csv",value=T)

	degs_in_mouse<-read.csv(file=input_file,sep=",",check.names=F)

	list_genes<-as.character(degs_in_mouse$Gene_name)

	statistics_degs<-degs_in_mouse[,c(1:7)]

	#matrix_mm<-degs_in_mouse[,-c(1:8)]
	
	#rownames(matrix_mm)<-paste(degs_in_mouse$ID,degs_in_mouse$Gene_name,sep=";")

	current_samples<- colnames(degs_in_mouse[,-c(1:7)])

	mm_sampleinfo_rid<-mm_sampleinfo[which(mm_sampleinfo$Sample_name%in%current_samples),]
	conds<-unique(mm_sampleinfo_rid$Condition2)
	filename<-as.character(mm_sampleinfo_rid$Sample_name)

	#
	# Read the matrix with the TPM values
	#
	#tpm_matrix<-fread("/mnt/data/lab/gmt_data/data_brain/adi_data_november/Adi_IDH_r123.tpm.csv",data.table=F)
	tpm_matrix<-degs_in_mouse[,-c(1:7)]
	#get all degs
	# tpm_matrix<-tpm_matrix[tpm_matrix$ID%in%degs_in_mouse$ID,]
	# #add rownames to the matrix
	rownames(tpm_matrix)<-paste(degs_in_mouse$ID,degs_in_mouse$Gene_name,sep=";")
	# #select the gene-expression profiles from samples in the current experiment
	matrix_mm2<-tpm_matrix[,which(colnames(tpm_matrix)%in%filename)]
	#order the matrix for the sampleinfo
	matrix_mm2<-matrix_mm2[,match(as.character(mm_sampleinfo_rid$Sample_name),colnames(matrix_mm2))]
	colnames(matrix_mm2)<-mm_sampleinfo_rid$Sample_name

	#matrix_mm2<-data.frame(matrix_mm[,match(mm_sampleinfo_rid$Sample_name,colnames(matrix_mm))])
	#colnames(matrix_mm2)<-colnames(matrix_mm)

	ann_mm_clean<-mm_sampleinfo_rid[,c("Condition","Condition2")]

	rownames(ann_mm_clean)<-mm_sampleinfo_rid$Sample_name

	set.seed(123) 
	
	#	
	# 4) Filter based on DEGs analysis 
	#

	#get the matrix with the gene-expression values
	statistics_degs<-degs_in_mouse[,c(1:8)]

	#select the differentially expressed genes
	selected_genes_with_degs<-as.character(statistics_degs[statistics_degs$padj<=pval & abs(statistics_degs$log2FoldChange)>=lfc,1])
	
	#select the up and down regulated genes
	up_genes<-statistics_degs[statistics_degs$padj<=pval & statistics_degs$log2FoldChange>=lfc,1]
	down_genes<-statistics_degs[statistics_degs$padj<=pval & statistics_degs$log2FoldChange <= -c(lfc),1]

	# save the DEGs in mm
	output_degs<-paste(paste(directory_exp,pval,lfc,sep="_"),".txt",sep="")

	output_txt<-paste(paste(output_dir,output_degs,sep="/"))
	
	write.table(selected_genes_with_degs,file=output_txt,sep="\t",row.names=T,quote=F)

	#obtain the ortologs genes with hs
	ortholog_table_mgi<-ortholog_table[which(ortholog_table[,5]%in%selected_genes_with_degs),]
	#
	# The definitive version of biomaRt is 102 -- the default, to change the version use the parameter "version"
	#
	print("WARNING:: The version of biomaRt used on 04/Feb/2021  (EU date) is 102, use convertMouseGeneList2 to use a old version ")
	print("WARNING:: The version of biomaRt used on 04/Feb/2021  (EU date) is 102, use convertMouseGeneList2 to use a old version")
	print("WARNING:: The version of biomaRt used on 04/Feb/2021  (EU date) is 102, use convertMouseGeneList2 to use a old version")
	print("WARNING:: The version of biomaRt used on 04/Feb/2021  (EU date) is 102, use convertMouseGeneList2 to use a old version")
	print("WARNING:: The version of biomaRt used on 04/Feb/2021  (EU date) is 102, use convertMouseGeneList2 to use a old version")
	
	ortholog_table_ensembl<-unique(convertMouseGeneList2(selected_genes_with_degs))

	ortholog_table_res<-integrated_conversion_mouse_human(mgi=ortholog_table_mgi,ensembl=ortholog_table_ensembl)
	print("WARNING:: I am using the ensembl version 102, 04/Feb/2021. We locked the analysis using this repository.")
	
	output_degs2<-paste(paste(directory_exp,pval,lfc,sep="_"),".ort.txt",sep="")
	output_txt2<-paste(paste(output_dir,output_degs2,sep="/"))
	
	#
	# The file for GSVA analysis is the list of ort. genes in human - this analysis works only for C1A and C1B not C1 and C2
	#
	
	file_gsva<-output_txt2
	
	hs_genes_sig<-unique(ortholog_table_res[,1])
	mm_genes_sig<-unique(ortholog_table_res[,2])

	hs_genes_sig_up<-ortholog_table_res[ortholog_table_res[,2]%in%up_genes,1]
 	hs_genes_sig_down<-ortholog_table_res[ortholog_table_res[,2]%in%down_genes,1]

	ortholog_table_res$status[ortholog_table_res[,1]%in%hs_genes_sig_up]<-"up"
	ortholog_table_res$status[ortholog_table_res[,1]%in%hs_genes_sig_down]<-"down"
	
	#export the ortologs genes 	
	write.table(ortholog_table_res,file=output_txt2,sep="\t",row.names=T,quote=F)

  #
  # Create matrix h.sapiens --NOT-- processed with combat
	#
	
	col_to_select<-c(colnames(input_ge[,c(1:6)]),"Type")
  matrix_to_plot_hs_nocombat<-t(input_ge[,-c(which(colnames(input_ge)%in%col_to_select))])
	colnames(matrix_to_plot_hs_nocombat)<-unlist(lapply(X=strsplit(input_ge$Sample,split="\\-"),FUN=function(X)paste(X[1:3],collapse="-")))

	matrix_to_plot_hs_nocombat2<-matrix_to_plot_hs_nocombat[which(rownames(matrix_to_plot_hs_nocombat)%in%hs_genes_sig),]

	#
	# Create matrix mouse with all genes derived from the DEGs analysis
	#
	        
	output_mm_boxplot<-paste("Mus_musculus",directory_exp,".method.fc",pval,lfc,".boxplot.allgenes.png",sep="")

	matrix_mm3<-data.frame(gs=unlist(lapply(strsplit(rownames(matrix_mm2),split=";"),"[[",2)),matrix_mm2)
	colnames(matrix_mm3)<-c("gs",colnames(matrix_mm2))

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
	
	###############################################
	# Plot human results
	################################################ 

	setwd(output_dir)

	annotation_all<-ann3[,c("Samples","type","IDH1","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]
	
	#
	# No Combat: all cancers
	#
	
	output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".nocombat.pdf",sep="")
  # output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".nocombat.var.pdf",sep="")

	plot_heatmap_hs(matrix_to_plot_hs_nocombat2,annotation_all,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=TRUE,number_clusters_row=2,png_output=F)	
        
	# list_degs<-split(ortholog_table_res,ortholog_table_res$status)
	# 
	# for(up_down_genes in 1:length(list_degs)){
	#  
	# current_degs_for_heatmap<-list_degs[[up_down_genes]][,1]
	# status_degs_for_heatmap<-names(list_degs)[[up_down_genes]]
	# 
	# output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",status_degs_for_heatmap,".nocombat.png",sep="")
	# output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",status_degs_for_heatmap,".nocombat.var.png",sep="")
	# if(length(current_degs_for_heatmap)>10){  
	# plot_heatmap_hs(matrix_to_plot_hs_nocombat2[which(rownames(matrix_to_plot_hs_nocombat2)%in%current_degs_for_heatmap),],annotation_all,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=TRUE,number_clusters_row=2)	
	# }
	# 
	# }
	
	#scores_hs<-CrossSpeciesComparison(mat_exp=matrix_to_plot_hs_nocombat2,ortholog_table_res,species="hs")
   	
    	#annotation_with_scores<-merge(annotation_all,scores_hs,by.x="Samples",by.y="ID")
   	
    	#output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".nocombat.withSCORE.pdf",sep="")
    	# output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".nocombat.var.withSCORE.png",sep="")
   	
    	#plot_heatmap_hs(matrix_to_plot_hs_nocombat2,annotation_with_scores,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=TRUE,number_clusters_row=2,png_output=F) 

	
	#
	# No Combat: only LGG samples
	#
	
 	annotation_all<-ann3[,c("Samples","type","IDH1","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]
	annotation_all_lgg<-annotation_all[annotation_all$type%in%"LGG",]
  	
	output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".LGG.nocombat.pdf",sep="")
	# output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".LGG.nocombat.var.png",sep="")
	
	matrix_to_plot_lgg<-matrix_to_plot_hs_nocombat2[,which(colnames(matrix_to_plot_hs_nocombat2)%in%annotation_all_lgg[,1])]

  plot_heatmap_hs(matrix_to_plot_lgg,annotation_all_lgg,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=TRUE,number_clusters_row=2,png_output=F)

  #
  # No Combat: only LGG samples up and down
  #
  
  # list_degs<-split(ortholog_table_res,ortholog_table_res$status)
  # 
  # for(up_down_genes in 1:length(list_degs)){
  #   
  #   current_degs_for_heatmap<-list_degs[[up_down_genes]][,1]
  #   status_degs_for_heatmap<-names(list_degs)[[up_down_genes]]
  #   
  #   output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",status_degs_for_heatmap,".LGG.nocombat.png",sep="")
  #   output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",status_degs_for_heatmap,".LGG.nocombat.var.png",sep="")
  # 
  # if(length(current_degs_for_heatmap)>10){  
  #     
  #   plot_heatmap_hs(matrix_to_plot_lgg[which(rownames(matrix_to_plot_lgg)%in%current_degs_for_heatmap),],annotation_all_lgg,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=TRUE,number_clusters_row=2)
  # }
  #   
  # }
  
	# Heatmap LGG with scores 

  #scores_hs_lgg<-CrossSpeciesComparison(mat_exp=matrix_to_plot_lgg,ortholog_table_res,species="hs")
  	       
 # annotation_with_scores<-merge(annotation_all,scores_hs_lgg,by.x="Samples",by.y="ID")
  
  #output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".LGG.nocombat.withSCORE.pdf",sep="")

  # output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".LGG.nocombat.var.withSCORE.png",sep="")
  			                
  #plot_heatmap_hs(matrix_to_plot_lgg,annotation_with_scores,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=TRUE,number_clusters_row=2,png_output=F) 
  
  # Heatmap samples IDH1 + TP53 + ATRX with scores 
  #annotation_all_compare_genes<-ann3[,c("Samples","type","IDH1","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]
  #annotation_all_compare_genes<-annotation_all_compare_genes[annotation_all_compare_genes$IDH1=="1" & annotation_all_compare_genes$type=="LGG",]
  
  #matrix_lgg_triple_mut<-matrix_to_plot_hs_nocombat2[,which(colnames(matrix_to_plot_hs_nocombat2)%in%annotation_all_compare_genes[,1])]
  
  #scores_lgg_triple_mut<-CrossSpeciesComparison(mat_exp=matrix_lgg_triple_mut,ortholog_table_res,species="hs")
  
  #annotation_lgg_triple_mut_with_scores<-merge(annotation_all_compare_genes,scores_lgg_triple_mut,by.x="Samples",by.y="ID")
  #annotation_lgg_triple_mut_with_scores<-annotation_lgg_triple_mut_with_scores[,-12]
  
  #output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".LGG_TRIPLEMUT.nocombat.withSCORE.pdf",sep="")
  
  # output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".LGG_TRIPLEMUT.nocombat.var.withSCORE.png",sep="")
  
  #plot_heatmap_hs(matrix_lgg_triple_mut,annotation_lgg_triple_mut_with_scores,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=TRUE,number_clusters_row=2,png_output=F) 
  
  #
  # Only for C1A and C1B plot the heatmap with TME IDH1 + TP53 + ATRX
  #
  
  merge_exp_with_TME = FALSE
  
  if(merge_exp_with_TME==TRUE){
    
  source("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/colors_scheme_chart2.R")
  
  load("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/df.tcga.groups.RData")
  load("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/TMEinfiltration_TCGA_GBMplusLGG.RData")
  
  annotation_all_compare_genes<-ann3[,c("Samples","type","IDH1","IDH1_R132H","IDH1_other","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]
  annotation_all_compare_genes$status2<-rep(NA,nrow(annotation_all_compare_genes))
  annotation_all_compare_genes$status2[-which((annotation_all_compare_genes$IDH1_R132H==1 | annotation_all_compare_genes$IDH1_other==1) & annotation_all_compare_genes$TP53==1 & annotation_all_compare_genes$ATRX==1)]<-0
  annotation_all_compare_genes$status2[which((annotation_all_compare_genes$IDH1_R132H==1 | annotation_all_compare_genes$IDH1_other==1) & annotation_all_compare_genes$TP53==1 & annotation_all_compare_genes$ATRX==1)]<-1
  annotation_all_compare_genes<-annotation_all_compare_genes[annotation_all_compare_genes$status2=="1" & annotation_all_compare_genes$type=="LGG",]
  
  df.tcga.groups[,1]<- unlist(lapply(strsplit(df.tcga.groups[,1],split="\\-"),FUN=function(X){paste(X[1:3],collapse="-")}))
  df.tcga.groups[,2]<-gsub(df.tcga.groups[,2],pattern="/",replacement="")
  
  annotation_all_compare_genes2<-merge(annotation_all_compare_genes,df.tcga.groups,by.x="Samples",by.y="Sample")
  
  groups_coxph<-read.delim(file="OS.Comp8_C1A_vs_C1B_0.05_1.gsva_MCP_0_LGG_TRIPLE_MUT.coxph_thrcox.0.01.up_dw.txt",sep="\t")
  annotation_all_compare_genes2<-merge(annotation_all_compare_genes2,groups_coxph[,c(1:2)],by.x="Samples",by.y="ID")
  
  matrix_lgg_triple_mut<-matrix_to_plot_hs_nocombat2[,which(colnames(matrix_to_plot_hs_nocombat2)%in%annotation_all_compare_genes2[,1])]
  
  scores_lgg_triple_mut<-CrossSpeciesComparison(mat_exp=matrix_lgg_triple_mut,ortholog_table_res,species="hs")
  
  annotation_lgg_triple_mut_with_scores<-merge(annotation_all_compare_genes2,scores_lgg_triple_mut,by.x="Samples",by.y="ID")
  annotation_lgg_triple_mut_with_scores<-annotation_lgg_triple_mut_with_scores[,-which(colnames(annotation_lgg_triple_mut_with_scores)%in%"status2")]
  
  output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".LGG_TRIPLEMUT.nocombat.withSCOREandTME.pdf",sep="")
  
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  
  mat_hs_cluster3<-log(matrix_lgg_triple_mut+1,2)
  mat_hs_cluster3ok<-t(apply(mat_hs_cluster3,1,cal_z_score))
  colnames(mat_hs_cluster3ok)<-colnames(mat_hs_cluster3)
  
  colAnn <- HeatmapAnnotation(df=annotation_lgg_triple_mut_with_scores[,-1], which="col",
                              col=colours2,
                              annotation_width=unit(c(1, 4), "cm"),
                              gap=unit(1, "mm"),
                              annotation_name_side = "left",
                              border=TRUE)
  
  output_integrateTME<-paste(paste("C1AC1B_integrate_TME",lfc,pval,sep="_"),".pdf",sep="")
  
  pdf(output_integrateTME,width=15,height=10)
  
  mat_hs_cluster3ok[is.na(mat_hs_cluster3ok)]<-0
  
  ht1<-Heatmap(mat_hs_cluster3ok,
               clustering_distance_rows = "euclidean",
               clustering_method_rows = "complete",
               clustering_distance_columns = "euclidean",
               clustering_method_columns = "complete",
               top_annotation=colAnn,
               column_km = 3,
               row_km = 2,
               show_row_names = FALSE,
               show_column_names = FALSE,
               height = 10
              )
  
  tme.TCGA2<-merge(annotation_lgg_triple_mut_with_scores,tme.TCGA,by.x="Samples",by.y="Case")
  tme.TCGA3<-tme.TCGA2[,-which(colnames(tme.TCGA2)%in%colnames(annotation_lgg_triple_mut_with_scores))]
  tme.TCGA3[,1]<-unlist(lapply(strsplit(tme.TCGA3[,1],split="\\-"),FUN=function(X){paste(X[1:3],collapse="-")}))
  
  tme_matrix<-t(tme.TCGA3[tme.TCGA3[,1]%in%colnames(mat_hs_cluster3ok),-c(1)])
  
  samplesID<-tme.TCGA3[tme.TCGA3[,1]%in%colnames(mat_hs_cluster3ok),1]
  colnames(tme_matrix)<-samplesID
  
  list_terms_to_use<-c("Bcells","Macrophages_M2","MDM_Consesus","Endothelial",
                       "Neutrophils","CD8_Tcells","NKcells","CD4_Tcells","T_regs",
                       "Microglia_Consensus","APC_activation","Macrophages_M1",
                       "T_cell_activation","Exhaustion","Memory_Tcells","MicrogliaMacrophagesCorrGenes")
  
  tme_matrix2<-tme_matrix[rownames(tme_matrix)%in%list_terms_to_use,]
  
  ht2<-Heatmap(tme_matrix2,
               clustering_distance_rows = "euclidean",
               clustering_method_rows = "complete",
               clustering_distance_columns = "euclidean",
               clustering_method_columns = "complete",
               show_row_names = TRUE,
               show_column_names = FALSE,
               height = 10,
               row_km = 2
  )
  
  htlist<-ht1 %v% ht2
  
  print(draw(htlist))
  
  #
  # Try to select in a different way the TME microenvironment terms
  #
  
  samples_in_km_column<-column_order(ht1)
  
  kDF<-data.frame()
  
  for(skcm in 1:length(samples_in_km_column)){
  
  current_cluster<-names(samples_in_km_column)[[skcm]]
  
  idx<-samples_in_km_column[[skcm]]
  samplesID<-colnames(mat_hs_cluster3ok)[idx]
  
  kDF<-rbind(kDF,data.frame(current_cluster,samplesID))
  
  }
  colnames(kDF)<-c("clusters_km_col","samples")
  
  kDFTCGAtme<-merge(kDF,tme.TCGA3,by.x="samples",by.y="Sample")
  
  ann_output_tme<-merge(kDF,annotation_lgg_triple_mut_with_scores,by.x="samples",by.y="Samples")
  output_integrateTME2<-paste(paste("C1AC1B_integrate_TME",lfc,pval,sep="_"),".txt",sep="")
  
  write.table(ann_output_tme,file=output_integrateTME2,sep="\t",row.names=F,quote=F)
  
  res_gain<-information.gain(clusters~., kDFTCGAtme[,-1])
  thr<-quantile(res_gain[,1])[2]
  idx_tme<-which(res_gain>thr)
  tme_to_save_FS<-rownames(res_gain)[idx_tme]
  
  tme_matrix2<-tme_matrix[rownames(tme_matrix)%in%tme_to_save_FS,]
  
  ht3<-Heatmap(tme_matrix2,
               clustering_distance_rows = "euclidean",
               clustering_method_rows = "complete",
               clustering_distance_columns = "euclidean",
               clustering_method_columns = "complete",
               show_row_names = TRUE,
               show_column_names = FALSE,
               height = 10
  )
  
  htlist2<-ht1 %v% ht3
  
  print(draw(htlist2))
  
  
  dev.off()
  
  }
  
	#
	# compare degs human with degs mouse
	#
	
	# get only the samples with mutations in IDH1_R132H, TP53, ATRX
	annotation_all_compare_genes<-ann3[,c("Samples","type","IDH1","IDH2","IDH1_other","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]
	annotation_all_compare_genes$status2<-rep(NA,nrow(annotation_all_compare_genes))
	annotation_all_compare_genes$status2[-which((annotation_all_compare_genes$IDH1==1 | annotation_all_compare_genes$IDH2==1) & annotation_all_compare_genes$TP53==1 & annotation_all_compare_genes$ATRX==1)]<-0
	annotation_all_compare_genes$status2[which((annotation_all_compare_genes$IDH1==1 | annotation_all_compare_genes$IDH2==1) & annotation_all_compare_genes$TP53==1 & annotation_all_compare_genes$ATRX==1)]<-1
	
	# matrix_to_compare_degs<-matrix_to_plot_hs_nocombat2[,which(colnames(matrix_to_plot_hs_nocombat2)%in%annotation_all_compare_genes[,1])]
	# 
	# matrix_to_compare_degs<-matrix_to_compare_degs[,match(annotation_all_compare_genes[,1],colnames(matrix_to_compare_degs))]
	
	mat_temp_trs<-t(matrix_to_plot_hs_nocombat2)
	mat_temp_trs_with_id<-data.frame(ID=rownames(mat_temp_trs),mat_temp_trs)
	
	size_ann<-ncol(annotation_all_compare_genes)
	
	merged_exp_ann<-merge(annotation_all_compare_genes,mat_temp_trs_with_id,by.x="Samples",by.y="ID")
	
	annotation_all_compare_genes2<-merged_exp_ann[,c(1:size_ann)]
	matrix_to_compare_degs<-t(merged_exp_ann[,-c(1:size_ann)])
	colnames(matrix_to_compare_degs)<-merged_exp_ann$Samples
	
	
	# if(directory_exp!="Comp8_C1A_vs_C1B"){
	#   
	# (matrix_hs=matrix_to_compare_degs,annotation_for_comp_degs=annotation_all_compare_genes,degs_mouse=ortholog_table_res,lfc=1,pvalue=0.05,output_string=paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".nocombat",sep=""))
	# 
	# }

	mlForModel(matrix_hs=matrix_to_compare_degs,
	           annotation_for_classes=annotation_all_compare_genes2,
	           matrix_all_genes=matrix_to_plot_hs_nocombat[-which(rownames(matrix_to_plot_hs_nocombat)%in%ortholog_table_res[,1]),],
	           output_string=paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".nocombat",sep=""),
	           randGenes=FALSE,
	           permutation=FALSE)

  # get only the samples with mutations in IDH1_R132H, TP53, ATRX only LGG
	annotation_all_compare_genes<-ann3[,c("Samples","type","IDH1","IDH2","IDH1_other","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]
	
	annotation_all_compare_genes$status2<-rep(NA,nrow(annotation_all_compare_genes))
	annotation_all_compare_genes$status2[-which((annotation_all_compare_genes$IDH1==1 | annotation_all_compare_genes$IDH2==1) & annotation_all_compare_genes$TP53==1 & annotation_all_compare_genes$ATRX==1)]<-0
	annotation_all_compare_genes$status2[which((annotation_all_compare_genes$IDH1==1 | annotation_all_compare_genes$IDH2==1) & annotation_all_compare_genes$TP53==1 & annotation_all_compare_genes$ATRX==1)]<-1
	
	annotation_all_compare_genes_lgg<-annotation_all_compare_genes[annotation_all_compare_genes$type=="LGG",]

	mat_temp_trs<-t(matrix_to_plot_hs_nocombat2)
	mat_temp_trs_with_id<-data.frame(ID=rownames(mat_temp_trs),mat_temp_trs)
	
	size_ann<-ncol(annotation_all_compare_genes)
	
	merged_exp_ann<-merge(annotation_all_compare_genes_lgg,mat_temp_trs_with_id,by.x="Samples",by.y="ID")
	
	annotation_all_compare_genes2<-merged_exp_ann[,c(1:size_ann)]
	matrix_to_compare_degs_lgg<-t(merged_exp_ann[,-c(1:size_ann)])
	colnames(matrix_to_compare_degs_lgg)<-merged_exp_ann$Samples
	
	# if(directory_exp!="Comp8_C1A_vs_C1B"){
	#   
	# compareDegsHs_MM(matrix_hs=matrix_to_compare_degs_lgg,annotation_for_comp_degs=annotation_all_compare_genes_lgg,degs_mouse=ortholog_table_res,lfc=1,pvalue=0.05,output_string=paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".nocombat_LGG",sep=""))
	# 
	# }

  mlForModel(matrix_hs=matrix_to_compare_degs_lgg,
             annotation_for_classes=annotation_all_compare_genes2,
             matrix_all_genes=matrix_to_plot_hs_nocombat[-which(rownames(matrix_to_plot_hs_nocombat)%in%ortholog_table_res[,1]),],
             output_string=paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".nocombat_LGG",sep=""),
             randGenes=FALSE,
             permutation=FALSE)


	
  # #
  # # get only the samples with mutation in IDH1 + TP53 + ATRX GBM + LGG
  # #
  # ann_TRIPLEMUT<-ann3[(ann3$IDH1==1) & (ann3$TP53==1 & ann3$ATRX==1),]
  # 
  # matrix_TRIPLEMUT<-matrix_to_plot_hs_nocombat2[,which(colnames(matrix_to_plot_hs_nocombat2)%in%ann3[,1])]
  
  
	#
	# Survival analysis
  #
  
  # define a function to average the expression of the genes in repeated samples
  average_samples<-function(X,samples_to_average){
    
    dfResAverage<-data.frame(1:nrow(X))
    
    for(asi in 1:length(samples_to_average)){
      
      current_asi<-samples_to_average[asi]
      
      Xrid<-X[,colnames(X)%in%current_asi]
      
      Xrid_avg<-data.frame(rowMeans(Xrid))
      
      dfResAverage<-cbind(dfResAverage,Xrid_avg)
      
    }
    
    dfResAverage2<-dfResAverage[,-1]
    colnames(dfResAverage2)<-samples_to_average
    
    Xwithoutrepeated_samples<-X[,-which(colnames(X)%in%samples_to_average)]
    
    dfFinal<-cbind(Xwithoutrepeated_samples,dfResAverage2)
    
    return(dfFinal)
  }
  
  # annotation_all<-ann3[,c("Samples","type","IDH1","IDH1_R132H","IDH1_other","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]
  #   
  #   if(directory_exp=="Comp8_C1A_vs_C1B"){
  #     
  # 	if(do_gsva==1){
  # 	
  # 	#get thhe samples with IDH1 + TP53 + ATRX mutations
  # 	annotation_mut<-annotation_all[annotation_all$type%in%"LGG" & (annotation_all$IDH1==1 & annotation_all$TP53==1 & annotation_all$ATRX==1),]
  # 	#get thhe samples with LGG and IDH1
  # 	annotation_only_idh<-annotation_all[annotation_all$type%in%"LGG" & (annotation_all$IDH1),]
  #   annotation_gbm_lgg<-annotation_all
  #   annotation_idh_codel_nocodel<-annotation_all[annotation_all$type%in%"LGG" & (annotation_all$IDH_codel_subtype=="IDHmut-codel" | annotation_all$IDH_codel_subtype=="IDHmut-non-codel"),]
  #   annotation_idh_codel<-annotation_all[annotation_all$type%in%"LGG" & (annotation_all$IDH_codel_subtype=="IDHmut-codel"),]
  #   annotation_idh_nocodel<-annotation_all[annotation_all$type%in%"LGG" & (annotation_all$IDH_codel_subtype=="IDHmut-non-codel"),]
  #   annotation_gbm_only<-annotation_all[annotation_all$type%in%"GBM",]
  # 	
  #   #
  #   # Define a function to prepare the input for coxph
  #   #
  #   
  #   prepareMatrixForCox<-function(mat_exp,annotation_file){
  #     
  #                         matrix_for_gsva<-log(mat_exp+1,2) #convert in log scale for gsva
  #                         
  #                         matrix_for_mut_temp<-matrix_for_gsva[,which(colnames(matrix_for_gsva)%in%annotation_file[,1])] #get the ge samples that are in the annotation file
  #                         nids<-data.frame(table(colnames(matrix_for_mut_temp))) #find the repeated samples
  #                         samples_to_average<-nids[nids[,2]>1,1] 
  #                         matrix_for_mut<-average_samples(matrix_for_mut_temp,samples_to_average) #collapse the expression values of the replicated samples
  #                         
  #                         return(matrix_for_mut)
  #                         
  #   }
  #   
  #   matrices_for_gsva<-vector(mode="list",7)
  #   names(matrices_for_gsva)<-c("LGG_TRIPLE_MUT","only_LGG_IDH","GBM_LGG","IDH_CODEL_NOCODEL","IDH_CODEL","IDH_NOCODEL","only_GBM")
  #   
  #   list_annotations<-c("annotation_mut","annotation_only_idh","annotation_gbm_lgg","annotation_idh_codel_nocodel","annotation_idh_codel","annotation_idh_nocodel","annotation_gbm_only")
  #   
  #   matrix_for_gsva<-matrix_to_plot_hs_nocombat
  #   
  #   for(lafm in 1:length(list_annotations)){
  #   
  #   annotation_temp_for_cox<-get(list_annotations[lafm])  
  #   mat_temp_for_cox<-prepareMatrixForCox(matrix_for_gsva,annotation_temp_for_cox)
  #   
  #   matrices_for_gsva[[lafm]]<-mat_temp_for_cox
  #   
  #   }
  #   
  #   
  #   
  #   
  #   # 
  #   # 	# get the matrix for gsva
  #   # 	matrix_for_gsva<-log(matrix_to_plot_hs_nocombat+1,2)
  #   # 	
  #   # 	matrix_for_mut_temp<-matrix_for_gsva[,which(colnames(matrix_for_gsva)%in%annotation_mut[,1])]
  #   # 	nids<-data.frame(table(colnames(matrix_for_mut_temp)))
  #   # 	samples_to_average<-nids[nids[,2]>1,1]
  #   # 	matrix_for_mut<-average_samples(matrix_for_mut_temp,samples_to_average)
  #   # 	
  #   # 	matrix_only_idh_temp<-matrix_for_gsva[,which(colnames(matrix_for_gsva)%in%annotation_only_idh[,1])]
  #   # 	nids<-data.frame(table(colnames(matrix_only_idh_temp)))
  #   # 	samples_to_average<-nids[nids[,2]>1,1]
  #   # 	matrix_only_idh<-average_samples(matrix_only_idh_temp,samples_to_average
  # 	# 
  # 	# matrices_for_gsva<-vector(mode="list",2)
  # 	# names(matrices_for_gsva)<-c("LGG_TRIPLE_MUT","only_LGG_IDH")
  # 	# #names(matrices_for_gsva)<-c("LGG_MUT","only_LGG_IDH1","GBM_LGG")
  # 	# 	
  # 	# matrices_for_gsva[[1]]<-matrix_for_mut
  # 	# matrices_for_gsva[[2]]<-matrix_only_idh
  # 	# #matrices_for_gsva[[3]]<-log(matrix_to_plot_hs_nocombat+1,2)
  #   
  #   #Download the clinical data for LGG and GBM
  #   clin.lgg <- GDCquery_clinic("TCGA-LGG", "clinical")
  #   clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
  #   
  #   #Download the important columns for the analysis
  #   features_to_select<-c("submitter_id","days_to_death","days_to_last_follow_up","vital_status")
  #   clin.lgg2<-clin.lgg[,colnames(clin.lgg)%in%features_to_select]
  #   clin.gbm2<-clin.gbm[,colnames(clin.gbm)%in%features_to_select]
  #   
  # 	for(an_gs in 1:length(matrices_for_gsva)){
  # 		
  # 	  # read the file with the differentially expressed genes           
  #     samplesc1ac1b<-gsub(sapply(strsplit(file_gsva,split="/"),"[[",13),pattern=".ort.txt",replacement="")
  #     
  #     # use the DEGs to compute a score with gsva
  # 		scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",max.sz=2000,cut_off_gsva=0)
  # 		
  # 		# extract the name of the experiment
  # 		current_score_gsva<-setDT(scores_gsva[[2]])
  # 	  	
  # 		output_string<-paste(scores_gsva[[1]],names(matrices_for_gsva)[[an_gs]],sep="_")
  # 		
  # 		output_pdf<- paste(paste("OS",sep="."),samplesc1ac1b,output_string,"pdf",sep=".")
  # 
  # 		output_os_txt<- paste(paste("OS",sep="."),samplesc1ac1b,output_string,"txt",sep=".")
  # 		write.table(scores_gsva,file=output_os_txt,sep="\t",row.names=F,quote=F)
  # 
  # 		pdf(output_pdf)
  # 
  # 		#
  # 		# Prepare clinical data
  # 		#
  # 
  #     #     clin.lgg <- GDCquery_clinic("TCGA-LGG", "clinical")
  #     #     clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
  #     # 
  #     # 		features_to_select<-c("submitter_id","days_to_death","days_to_last_follow_up","vital_status")
  #     #     clin.lgg2<-clin.lgg[,colnames(clin.lgg)%in%features_to_select]
  #     #     clin.gbm2<-clin.gbm[,colnames(clin.gbm)%in%features_to_select]
  #     
  # 		# put toghether the clinical data
  #     input_os_gsva_all<-rbind(clin.lgg2,clin.gbm2)
  #     
  # 		# input_os_gsva_lgg<-merge(clin.lgg2,current_score_gsva,by.x="submitter_id",by.y="ID")
  #     
  #     # use only the clinical data selected available in the matrix with the scores
  #     input_os_gsva_curr_group<-merge(input_os_gsva_all,current_score_gsva,by.x="submitter_id",by.y="ID")
  #     
  # 		#
  # 		# os with cut-off approach 
  # 		#
  # 				                
  # 		# plotSurvival(TCGA_surv=input_os_gsva_lgg,status="status_C1AC1B",title=paste(output_string,"custom thr. up_and_down"))
  #     plotSurvival(TCGA_surv=input_os_gsva_curr_group,status="status_C1AC1B",title=paste(output_string,"custom thr. up_and_down"))
  # 		
  # 		#
  # 		# os with best cut-point
  # 		#
  # 		
  # 		                
  #  		# res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="global_score",suffix_output=paste(output_string,"up_and_down"))
  #     res_os<-plotSurvival2(TCGA_surv=input_os_gsva_curr_group,score_to_use="global_score",suffix_output=paste(output_string,"up_and_down"))
  #     
  #  		#res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="down_score",suffix_output=paste(output_string,"only_down"))
  # 		                		
  #  		#res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="up_score",suffix_output=paste(output_string,"only_up"))
  # 
  # 		dev.off()
  # 	}
  # 
  # 	}
  #         
  #   #
  #   # Survival analysis coxph 
  #   #
  #   system("rm all_results_coxph_adi.txt")
  #   
  #   clin.lgg <- GDCquery_clinic("TCGA-LGG", "clinical")
  #   clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
  #   
  #   clin.lgg2<-clin.lgg[,colnames(clin.lgg)%in%features_to_select]
  #   clin.gbm2<-clin.gbm[,colnames(clin.gbm)%in%features_to_select]
  #   
  #   for(an_gs in 1:length(matrices_for_gsva[[1]])){
  #           
  #           features_to_select<-c("submitter_id","days_to_death","days_to_last_follow_up","vital_status")
  #           #this clinical data contains 515 patients, but they are merged with the selected samples in the gene-expression (function coxphAnalysis)
  #           # TCGA_surv_coxph<-clin.lgg2[,colnames(clin.lgg2)%in%features_to_select]
  #           TCGA_surv_coxph_all<-rbind(clin.lgg2,clin.gbm2)
  #           
  #           # get the annotation data derived previosuly to analyze
  #           temp_ann_coxph<-get(list_annotations[an_gs])
  #           
  #           #use only the samples in the annotation file 
  #           TCGA_surv_coxph<-TCGA_surv_coxph_all[which(TCGA_surv_coxph_all[,1]%in%temp_ann_coxph[,1]),]
  #           
  #           #define threhsolds  for the logRank analysis
  #           list_thr_cox<-c(0.01,0.001,0.05)
  # 
  #           for(thr_cox in list_thr_cox){
  #           
  #           # find the important genes with Coxph
  #           resCox<-coxphAnalysis(TCGA_surv_coxph,matrices_for_gsva[[an_gs]],file_gsva,thr_cox=thr_cox)
  #           resCox_genes<-resCox$Variable
  #           
  #           mtx_for_gsva_coxgenes<-matrices_for_gsva[[an_gs]]
  #           mtx_for_gsva_coxgenes2<-mtx_for_gsva_coxgenes[rownames(mtx_for_gsva_coxgenes)%in%resCox_genes,]
  #           
  #           samplesc1ac1b<-gsub(sapply(strsplit(file_gsva,split="/"),"[[",13),pattern=".ort.txt",replacement="")
  #           
  #           # define the type of gene list to use
  #           list_degs_cox<-c("up_dw","up","dw")
  # 
  #           all_res_cox2<-data.frame()
  #           
  #           for(type_degs_cox in list_degs_cox){
  #           
  #           print(paste("Analysis:",names(matrices_for_gsva)[[an_gs]],"-",thr_cox,"-",type_degs_cox))
  #             
  #           if(type_degs_cox == "up_dw"){
  #           
  # 	        print(type_degs_cox)
  # 
  #           scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes,cut_off_gsva=0)
  #           
  #           } 
  #             
  #           if(type_degs_cox == "up"){
  #           
  # 	        print(type_degs_cox)
  # 
  #           # select only up-regulated genes 
  #           resCox_genes_up<-as.character(resCox[resCox$status %in% "up",]$Variable)
  #           
  #           #if there are not sufficient up regulated genes from the cox analysis use directly the up regulated genes from the DEGs analysis
  #           if(length(resCox_genes_up)==0){
  #           
  #           resCox_genes_up<-as.character(ortholog_table_res[ortholog_table_res[,3]%in%"up",1])
  #             
  #           scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_up,cut_off_gsva=0)
  #               
  #           } else{
  #             
  #           # use the up regulated genes from coxph   
  #           resCox_genes_up<-as.character(resCox[resCox$status %in% "up",]$Variable)
  #             
  #           scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_up,cut_off_gsva=0)
  #             
  #           }
  #           
  #           }
  #             
  #           if(type_degs_cox == "dw"){
  # 
  # 	        print(type_degs_cox)
  # 
  #           # select only down-regulated genes
  #           resCox_genes_dw<-as.character(resCox[resCox$status %in% "dw",]$Variable)
  #           
  #           #if there are not sufficient down regulated genes from the cox analysis use directly the down regulated genes from the DEGs analysis
  #           if(length(resCox_genes_dw)==0){
  #             
  #           resCox_genes_dw<-as.character(ortholog_table_res[ortholog_table_res[,3]%in%"down",1])
  #             
  #           scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_dw,cut_off_gsva=0)
  #           
  #           } else {
  #             
  #           # use the down regulated genes from coxph   
  #           resCox_genes_dw<-as.character(resCox[resCox$status %in% "down",]$Variable)
  #           
  #           scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_dw,cut_off_gsva=0)
  #             
  #           } 
  #           
  #           } 
  #             
  #           current_score_gsva<-scores_gsva[[2]]
  #           
  #           output_string<-paste(scores_gsva[[1]],names(matrices_for_gsva)[[an_gs]],sep="_")
  #           
  #           output_pdf<- paste(paste("OS",sep="."),samplesc1ac1b,output_string,"coxph_thrcox",thr_cox,type_degs_cox,"pdf",sep=".")
  #           
  # 
  #           output_cox<-paste(paste("COX_RESULTS_thr_cox",thr_cox,sep="."),samplesc1ac1b,output_string,"coxph.txt",sep=".")
  #           write.table(resCox,file=output_cox,sep="\t",row.names=F,quote=F)
  # 
  #           output_os_txt<- paste(paste("OS",sep="."),samplesc1ac1b,output_string,"coxph_thrcox",thr_cox,type_degs_cox,"txt",sep=".")
  #           write.table(current_score_gsva,file=output_os_txt,sep="\t",row.names=F,quote=F)
  #           
  #           outpdf2<-paste(getwd(),output_pdf,sep="/")
  #           outtxt2<-paste(getwd(),output_cox,sep="/")
  #           
  #           line=paste(names(matrices_for_gsva)[[an_gs]],outpdf2,outtxt2,thr_cox,type_degs_cox,sep="\t")
  #           
  #           write(line,file="all_results_coxph_adi.txt",append=TRUE)
  #           
  #           pdf(output_pdf)
  #           
  #           #
  #           # Prepare clinical data
  #           #
  #           
  #           # input_os_gsva_lgg<-merge(clin.lgg2,current_score_gsva,by.x="submitter_id",by.y="ID")
  #           input_os_from_cox<-merge(TCGA_surv_coxph,current_score_gsva,by.x="submitter_id",by.y="ID",all.x=F,all.y=T)
  #           
  #           #
  #           # os with cut-off approach - up and down genes 
  #           #
  #           
  #           # plotSurvival(TCGA_surv=input_os_gsva_lgg,status="status_C1AC1B",title=paste(output_string,"custom thr. up_and_down"))
  #           plotSurvival(TCGA_surv=input_os_from_cox,status="status_C1AC1B",title=paste(output_string,"custom thr",type_degs_cox))
  #           
  #           #
  #           # os with best cut-point - up and down
  #           #
  #           
  #           # res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="global_score",suffix_output=paste(output_string,"up_and_down"))
  #           res_os<-plotSurvival2(TCGA_surv=input_os_from_cox,score_to_use="global_score",suffix_output=paste(output_string,type_degs_cox))
  #           
  #           #res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="down_score",suffix_output=paste(output_string,"only_down"))
  #           
  #           #res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="up_score",suffix_output=paste(output_string,"only_up"))
  #           
  #           dev.off()
  #           
  #           }
  #           
  #           }
  #           
  #         }
  #         
  #   }
  
  #
  # Create reports survival analysis
  #
 
  output_file<-paste("report_survival_analysis_coxph_logfc=",lfc,"qvalue=",pval,"html",sep=".")

  #rmarkdown::render('/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/res_04_02_2021/template_OS_Adi.Rmd', output_file = "report_survival_analysis_coxph_02_24_2021.html")

  #rmarkdown::render('/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/res_04_02_2021/template_OS_Adi.Rmd', output_file = output_file)

  # ALL_GROUPS<-NULL
  # ALL_PDF<-NULL
  # ALL_TXT<-NULL
  # ALL_PVALUE<-NULL
  # ALL_STATUS<-NULL
  # 
  # for(an_gs in 1:length(matrices_for_gsva)){
  #   
  #   list_thr_cox<-c(0.01,0.001,0.05)
  #   
  #   for(thr_cox in list_thr_cox){
  #     
  #     for(type_degs_cox in list_degs_cox){
  #       
  #       
  #       current_score_gsva<-scores_gsva[[2]]
  #       
  #       output_string<-paste(scores_gsva[[1]],names(matrices_for_gsva)[[an_gs]],sep="_")
  #       
  #       # output_pdf<- paste(getwd(),paste(paste("OS",sep="."),samplesc1ac1b,names(matrices_for_gsva)[[an_gs]],"coxph_thrcox",thr_cox,type_degs_cox,"pdf",sep="."),sep="/")
  #       output_pdf<- paste(getwd(),paste("OS",sep="."),samplesc1ac1b,output_string,"coxph_thrcox",thr_cox,type_degs_cox,"pdf",sep=".")
  #       
  #       output_cox<-paste(paste("COX_RESULTS_thr_cox",thr_cox,sep="."),samplesc1ac1b,output_string,"coxph.txt",sep=".")
  #       
  #       ALL_GROUPS<-c(ALL_GROUPS,names(matrices_for_gsva)[[an_gs]])
  #       ALL_PDF<-c(ALL_PDF,output_pdf)
  #       ALL_TXT<-c(ALL_TXT,output_pdf)
  #       ALL_PVALUE<-c(ALL_PVALUE,thr_cox)
  #       ALL_STATUS<-c(ALL_STATUS,type_degs_cox)
  #       
  #     }
  #   }
  #   
  #   
  # }
  # 
  # dfCoxphALL<-data.frame(ALL_GROUPS,
  #                        ALL_PDF,
  #                        ALL_TXT,
  #                        ALL_PVALUE,
  #                        ALL_STATUS)
  # 
  # write.table(dfCoxphALL,file="all_results_coxph_adi.txt",sep="\t",row.name=F,quote=F)
  
	#                       
	# Export filtered genes
  #                       
	        
	# sample_for_analysis<-annotation_all[,1]
	# 
	# matrix_to_plot_hs<-matrix_to_plot_hs_nocombat2
	# 
	# matrix_to_plot_hs<-matrix_to_plot_hs[,which(colnames(matrix_to_plot_hs)%in%sample_for_analysis)]
	# 
	# matrix_to_plot_hs<-matrix_to_plot_hs[,match(colnames(matrix_to_plot_hs),annotation_all[,1])]
  
  # 	mat_temp_trs<-t(matrix_to_plot_hs_nocombat2)
  # 	mat_temp_trs_with_id<-data.frame(ID=colnames(matrix_to_plot_hs_nocombat2),mat_temp_trs)
  # 	
  # 	size_ann<-ncol(annotation_all)
  # 	
  # 	merged_exp_ann<-merge(annotation_all,mat_temp_trs_with_id,by.x="Samples",by.y="ID")
  # 	
  # 	annotation_all<-merged_exp_ann[,c(1:size_ann)]
  # 	matrix_to_plot_hs<-t(merged_exp_ann[,-c(1:size_ann)])
  # 	
  # 	variance_genes<-rowVars(matrix_to_plot_hs)
  # 	variance_genes_sorted<-variance_genes[order(variance_genes)]
  # 	cut_off<-quantile(variance_genes_sorted)[2]
  # 
  # 	genes_after_variance_filter<-names(variance_genes_sorted[variance_genes_sorted>=cut_off])
  # 
  #   ortholog_table_res2<-ortholog_table_res[which(ortholog_table_res$hs%in%genes_after_variance_filter),]
  #  	output_degs2<-paste(paste(directory_exp,pval,lfc,sep="_"),".ort.var.txt",sep="")
  #   output_txt2_var<-paste(paste(output_dir,output_degs2,sep="/"))
  # 
  # 	write.table(ortholog_table_res2,file=output_txt2_var,sep="\t",row.names=T,quote=F)

	#
	# Feature selection by clusters
	#
	
	if(feature_selection == 1){

	cal_z_score <- function(x){
		      (x - mean(x)) / sd(x)
	}

	matrix_to_fs<-log(matrix_to_plot_hs_nocombat2+1,2)
	matrix_to_fs<-apply(matrix_to_fs,1,cal_z_score)
	matrix_to_fs[is.na(matrix_to_fs)]<-0
	
	library(mclust)
	library(FSelector)

	set.seed(123)
	mod1 <- Mclust(matrix_to_fs,mkcluster)
	mod2_classification<- mod1$classification
	
	matrix_to_fs2<-data.frame(clusters=mod2_classification, matrix_to_fs)

	res_gain<-information.gain(clusters~., matrix_to_fs2)
	thr<-quantile(res_gain[,1])[3]
	idx_genes<-which(res_gain>thr)
	genes_to_save_FS<-rownames(res_gain)[idx_genes]

	# GBM + LGG
	matrix_to_plot_fs<-matrix_to_plot_hs_nocombat2[rownames(matrix_to_plot_hs_nocombat2)%in%genes_to_save_FS,]

	output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS.pdf",sep="")
	# output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.var.FS.png",sep="")
  
	# check_ann<-merge(annotation_all,mod2_classification,by.x="Samples",by.y="")
	
	# output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.CHECKANN.FS.txt",sep="")
	# write.table(check_ann,file=output_hs,sep="\t",row.names=F,quote=F)
	mod_df<-data.frame(ID=rownames(matrix_to_fs),mod2_classification)
	mod_df$mod2_classification<-paste("cl",mod_df$mod2_classification,sep="")
  annotation_all2<-merge(annotation_all[,c("Samples","type","IDH1","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")],mod_df,by.x="Samples",by.y="ID")
    
	plot_heatmap_hs(matrix_to_plot_fs,annotation_all2,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=7,height=7,split=TRUE,km=FALSE,number_clusters_row=3,how_split="mod2_classification",strong_cluster=T,png_output=F)

	matrix_to_plot_fs<-matrix_to_plot_hs_nocombat2[rownames(matrix_to_plot_hs_nocombat2)%in%genes_to_save_FS,]
	
	# 
	# Heatmap for manuscript
	# 
	
	output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS_version2.pdf",sep="")
	# output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.var.FS_version2.pdf",sep="")
	
	# annotation_all2$mod2_classification<-paste("cl",annotation_all2$mod2_classification,sep="")
  annotation_all2$mod2_classification<-factor(as.character(annotation_all2$mod2_classification),levels=c("cl1","cl4","cl3","cl2"))

	plot_heatmap_hs(matrix_to_plot_fs,annotation_all2,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=7,height=7,split=TRUE,km=FALSE,number_clusters_row=3,how_split="mod2_classification",strong_cluster=T,png_output=F)

	# #
	# # GBM + LGG check no split
	# #
	# output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS.check.png",sep="")
	# 
	# # output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.var.FS.check.png",sep="")
	# 
	# ann_check<-cbind(annotation_all,mod2_classification)
	# 
	# table(ann_check$IDH1_R132H,ann_check$mod2_classification)
	# table(ann_check$clusters,ann_check$mod2_classification)
	# 
	# plot_heatmap_hs(matrix_to_plot_fs,ann_check,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=FALSE,number_clusters_row=3,strong_cluster=T,png_output=F)
	
  #   #
  # 	# Check cluster 3
  # 	#
  # 
  # 	output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS.check.CLUSTER3.png",sep="")
  # 
  # 	output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.var.FS.check.CLUSTER3.png",sep="")
  # 
  # 	ann_temp<-cbind(annotation_all,mod2_classification)
  # 	ann_check_cluster3<-unique(ann_temp[which(ann_temp$mod2_classification %in% 3),])
  # 	
  # 	# this sample is repeated two time i take the row in which the annotation is better
  # 	#idx2<- which(ann_check_cluster3[,1]%in%"TCGA-TQ-A7RK")[2]
  # 
  # 	#ann_check_cluster3<-ann_check_cluster3[-idx2,]
  # 	matrix_to_plot_fs_cluster3<-matrix_to_plot_fs[,colnames(matrix_to_plot_fs)%in%ann_check_cluster3[,1]]
  # 	temp_matrix_cluster3<-data.frame(samples=rownames(t(matrix_to_plot_fs_cluster3)),t(matrix_to_plot_fs_cluster3))
  # 	mat_exp_cluster3_withann<-merge(ann_check_cluster3,temp_matrix_cluster3,by.x="Samples",by.y="samples")
  # 
  # 	ann_cluster3_final<-mat_exp_cluster3_withann[,c(1:13)]
  # 	mat_exp_cluster3_final<-t(mat_exp_cluster3_withann[,-c(1:14)])	
  # 	colnames(mat_exp_cluster3_final)<-mat_exp_cluster3_withann[,1]
  # 
  # 	#matrix_to_plot_fs_cluster3<-matrix_to_plot_fs[,which(colnames(matrix_to_plot_fs)%in%ann_check_cluster3[,1])]
  # 
  # 	plot_heatmap_hs(mat_exp_cluster3_final,ann_cluster3_final,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=FALSE,number_clusters_row=3,how_split=ann_cluster3_final$ATRX,strong_cluster=T)
# 
# 	source("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/colors_scheme_chart.R")
# 
# 	pdf("manual_code_check_cluster3.pdf",width=15)
# 	mat_hs_cluster3<-log(mat_exp_cluster3_final+1,2)
# 	mat_hs_cluster3ok<-t(apply(mat_hs_cluster3,1,cal_z_score))
# 	colnames(mat_hs_cluster3ok)<-colnames(mat_hs_cluster3)
# 	
# 	colAnn <- HeatmapAnnotation(df=ann_cluster3_final[,-1], which="col",
# 				    col=colours,
# 				    annotation_width=unit(c(1, 4), "cm"),
# 				    gap=unit(1, "mm"),
# 				    annotation_name_side = "left",
# 				    border=TRUE)
# 
# 	ht1<-Heatmap(mat_hs_cluster3ok,
# 		clustering_distance_rows = "euclidean",
# 		clustering_method_rows = "complete",
# 		clustering_distance_columns = "euclidean",
# 		clustering_method_columns = "complete",
# 		top_annotation=colAnn,
# 		column_split=ann_cluster3_final$ATRX)
# 
# 	ht1ok<-draw(ht1,column_title="Stratification ATRX")
# 
#         ht2<-Heatmap(mat_hs_cluster3ok,
# 		     clustering_distance_rows = "euclidean",				
# 		     clustering_method_rows = "complete",							    
#    		     clustering_distance_columns = "euclidean",							
# 		     clustering_method_columns = "complete",							
# 		     top_annotation=colAnn)
# 
# 	ht2ok<-draw(ht2,column_title="NO Stratification ATRX")
# 
# 	rm(colours)
# 	dev.off()

# GBM + LGG Add scores C1 / C2
# 	scores_hs_fs<-CrossSpeciesComparison(mat_exp=matrix_to_plot_fs,ortholog_table_res,species="hs")
# 	              	
# 	annotation_with_scores<-cbind(annotation_all,scores_hs_fs[,-1])
# 		        
# 	output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".nocombat.FS.withSCORE.png",sep="")
# 			              	
# 	output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".nocombat.var.FS.withSCORE.png",sep="")
# 
#   plot_heatmap_hs(matrix_to_plot_fs,annotation_with_scores,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=FALSE,number_clusters_row=3,how_split=mod2_classification,strong_cluster=T)
        
	# LGG
  # 
  # 	annotation_fs_lgg<-annotation_all2[annotation_all2$type%in%"LGG",]
  # 	output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS_LGG.pdf",sep="")
  # 	# output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.var.FS_LGG.pdf",sep="")
  # 
  #   matrix_to_plot_fs_lgg<-matrix_to_plot_fs[,which(colnames(matrix_to_plot_fs)%in%annotation_fs_lgg[,1])]
  #   
  # 	plot_heatmap_hs(matrix_to_plot_fs_lgg,annotation_fs_lgg,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in LGG,no combat",width=10,split=TRUE,km=FALSE,number_clusters_row=3,how_split="mod2_classification",strong_cluster=T,png_output=F)
	
	# 
	# Heatmap for manuscript
	# 
	
	annotation_fs_lgg<-annotation_all2[which(annotation_all2$type%in%"LGG"),]
	
	matrix_to_plot_fs_lgg<-matrix_to_plot_fs[,which(colnames(matrix_to_plot_fs)%in%annotation_fs_lgg[,1])]
	
	# annotation_fs_lgg$mod2_classification<-factor(annotation_fs_lgg$mod2_classification,levels=c(1,3,4,2))
	
	output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.CHECKANN.FS.txt",sep="")
	write.table(annotation_fs_lgg,file=output_hs,sep="\t",row.names=F,quote=F)
	
	output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS_LGG_version2.pdf",sep="")
	plot_heatmap_hs(matrix_to_plot_fs_lgg,annotation_fs_lgg,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in LGG,no combat",width=8,height=8,split=TRUE,km=FALSE,number_clusters_row=3,how_split="mod2_classification",strong_cluster=T,png_output=F)
	
	# LGG Add scores C1 / C2
		
	# scores_hs_fs_lgg<-CrossSpeciesComparison(mat_exp=matrix_to_plot_fs_lgg,ortholog_table_res,species="hs")	
	# annotation_with_scores_fs_lgg<-cbind(annotation_fs_lgg,scores_hs_fs_lgg[,-1])
	# 
	# output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS_LGG.withSCORE.png",sep="")	
	# output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.var.FS_LGG.withSCORE.png",sep="")
	# 
	# matrix_to_plot_fs_lgg<-matrix_to_plot_fs[,which(colnames(matrix_to_plot_fs)%in%annotation_fs_lgg[,1])]
	#        
	# plot_heatmap_hs(matrix_to_plot_fs_lgg,annotation_with_scores_fs_lgg,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in LGG,no combat",width=10,split=TRUE,km=FALSE,number_clusters_row=3,how_split=annotation_fs_lgg$mod2_classification,strong_cluster=T)
	# 
	#
	#LGG Add scores C1/C2 FS plot by CODEL 
	#

	# output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS_LGG.withSCORE.CODEL.png",sep="")     	
	# output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.var.FS_LGG.withSCORE.CODEL.png",sep="")
	# 
	# plot_heatmap_hs(matrix_to_plot_fs_lgg,annotation_with_scores_fs_lgg,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in LGG,no combat",width=10,split=TRUE,km=FALSE,number_clusters_row=2,how_split=annotation_fs_lgg$IDH_codel_subtype,strong_cluster=T)

	# ml on IDH1 TP53 ATRX
	mat_temp_trs<-t(matrix_to_plot_fs)
	mat_temp_trs_with_id<-data.frame(ID=rownames(mat_temp_trs),mat_temp_trs)
	
	ann_temp<-ann3[,c("Samples","type","IDH1","IDH2","IDH1_R132H","IDH1_other","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]
	ann_temp$status2<-rep(NA,nrow(ann_temp))
	ann_temp$status2[-which((ann_temp$IDH1==1 | ann_temp$IDH2==1) & ann_temp$TP53==1 & ann_temp$ATRX==1)]<-0
	ann_temp$status2[which((ann_temp$IDH1==1 | ann_temp$IDH2==1) & ann_temp$TP53==1 & ann_temp$ATRX==1)]<-1
	
	size_ann<-ncol(ann_temp)
	
	merged_exp_ann<-merge(ann_temp,mat_temp_trs_with_id,by.x="Samples",by.y="ID")
	
	annotation_ml<-merged_exp_ann[,c(1:size_ann)]
	input_ml_fs<-t(merged_exp_ann[,-c(1:size_ann)])
	
	if(directory_exp!="Comp8_C1A_vs_C1B"){
	  
  	mlForModel(matrix_hs=input_ml_fs,
             annotation_for_classes=annotation_ml,
             matrix_all_genes=matrix_to_plot_hs_nocombat[-which(rownames(matrix_to_plot_hs_nocombat)%in%ortholog_table_res[,1]),],#remove informative genes from the entire matrix
             output_string=paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS",sep=""),
             randGenes=FALSE,
             permutation=FALSE)
	} 
	
	# ml on IDH1 TP53 ATRX only LGG        
	mat_temp_trs<-t(matrix_to_plot_fs)
	mat_temp_trs_with_id<-data.frame(ID=rownames(mat_temp_trs),mat_temp_trs)
	
	ann_temp<-ann3[,c("Samples","type","IDH1","IDH2","IDH1_R132H","IDH1_other","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]
	ann_temp$status2<-rep(NA,nrow(ann_temp))
	ann_temp$status2[-which((ann_temp$IDH1==1 | ann_temp$IDH2==1) & ann_temp$TP53==1 & ann_temp$ATRX==1)]<-0
	ann_temp$status2[which((ann_temp$IDH1==1 | ann_temp$IDH2==1) & ann_temp$TP53==1 & ann_temp$ATRX==1)]<-1
	
	ann_temp<-ann_temp[which(ann_temp$type%in%"LGG"),]
	
	size_ann<-ncol(ann_temp)
	
	merged_exp_ann<-merge(ann_temp,mat_temp_trs_with_id,by.x="Samples",by.y="ID")
	
	annotation_ml_lgg<-merged_exp_ann[,c(1:size_ann)]
	input_ml_fs<-t(merged_exp_ann[,-c(1:size_ann)])
	
	if(directory_exp!="Comp8_C1A_vs_C1B"){
	  
	mlForModel(matrix_hs=input_ml_fs,
	           annotation_for_classes=annotation_ml_lgg,
	           matrix_all_genes=matrix_to_plot_hs_nocombat[-which(rownames(matrix_to_plot_hs_nocombat)%in%ortholog_table_res[,1]),],
	           output_string=paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS_LGG",sep=""),
	           randGenes=FALSE,
	           permutation=FALSE,
	           iter_perm=10000)
	
	}
	
	ann_r132h_lgg<-ann3[,c("Samples","type","IDH1_R132H","IDH1_other","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]
	ann_r132h_lgg<-ann_r132h_lgg[ann_r132h_lgg$type=="LGG" & (ann_r132h_lgg$IDH1_R132H==1 | ann_r132h_lgg$IDH1_other==1),]
	
	matrix_onlyR132H_lgg<-matrix_to_plot_hs_nocombat2[,which(colnames(matrix_to_plot_hs_nocombat2)%in%ann_r132h_lgg[,1])]
	
	#
  # 	library(TCGAbiolinks)
  # 	clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
  # 	clin.lgg <- GDCquery_clinic("TCGA-LGG", "clinical")
  # 	
  # 	features_to_select<-c("submitter_id","days_to_death","days_to_last_follow_up","vital_status")
  # 	clin.gbm2<-clin.gbm[,colnames(clin.gbm)%in%features_to_select]
  # 	clin.lgg2<-clin.lgg[,colnames(clin.lgg)%in%features_to_select]
  # 
  # 	clin.tot<-rbind(clin.gbm2,clin.lgg2)
  # 	
  # 	mclust_classification<-data.frame(as.character(names(mod2_classification)),data.frame(mod2_classification))
  # 	colnames(mclust_classification)<-c("samples","mclust_clusters")
  # 
  # 	clin.tot2<-merge(clin.tot,mclust_classification,by.x="submitter_id",by.y="samples")
  # 
  #   	output_hs_os<-paste("Overall_survival",".",paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS.pdf",sep="")
  # 
  # 	# OS FS + LGG + GBM
  # 	pdf(output_hs_os)
  # 	plotSurvival(TCGA_surv=clin.tot2,status="mclust_clusters",title="OS GBM + LGG,clusters")
  # 	dev.off()
  # 
  # 	# OS FS + LGG
  #   output_hs_os<-paste("Overall_survival_LGG",".",paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS.pdf",sep="")
  # 
  # 	pdf(output_hs_os)
  # 	lgg_samples<-annotation_all[annotation_all$type=="LGG",1]
  # 	clin.tot2_lgg<-clin.tot2[clin.tot2[,1]%in%lgg_samples,]
  #   	plotSurvival(TCGA_surv=clin.tot2_lgg,status="mclust_clusters",title="OS LGG,clusters")
  # 	dev.off()

	# # OS FS + GSVA + LGG
	# genes_for_gsva<-read.delim(file_gsva)[,1]
	# ovp_gsva_fs<-length(intersect(genes_for_gsva,rownames(matrix_to_plot_fs)))

	# print(paste("The overlap between C1AC1B and FS is:",ovp_gsva_fs))
	
	}

  # 	#
  # 	# FS + TME 
  # 	#
  # 	if(feature_selection==1){
  # 
  # 	annotation_tme<-annotation_all
  # 
  # 	load("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_17_08_2020/msecrier_TME_adi/df.tcga.groups.RData")
  # 	colnames(df.tcga.groups)[c(1:2)]<-c("Samples","TME")
  # 	df.tcga.groups[,1]<- unlist(lapply(strsplit(df.tcga.groups[,1],split="\\-"),FUN=function(X){paste(X[1:3],collapse="-")}))
  # 
  # 	ann_tme<-merge(annotation_all,df.tcga.groups,by.x="Samples",by.y="Samples")
  # 
  # 	idx<-which(colnames(matrix_to_plot_fs) %in% ann_tme[,1])
  # 
  # 	matrix_to_plot_fs_TME<-matrix_to_plot_fs[,idx]
  # 
  #   	output_hs<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.FS.TME.png",sep="")
  #   	output_hs_var<-paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".nocombat.var.FS.TME.png",sep="")
  # 
  # 	mod2_classification_tme<-mod2_classification[which(names(mod2_classification)%in%ann_tme[,1])]
  # 
  # 	plot_heatmap_hs(matrix_to_plot_fs_TME,ann_tme,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=FALSE,number_clusters_row=3,how_split=mod2_classification_tme,strong_cluster=T)
  # 	
  # 	}

	#
	# Before to run the analysis with the methylation data, I export the gene-lists as xlsx
	#

 	output_degs<-paste(paste("REPORT",directory_exp,pval,lfc,sep="_"),"nocombat","nclust",mkcluster,"xlsx",sep=".")
  output_xlsx<-paste(paste(output_dir,output_degs,sep="/"))

	# genes DEGs m.musculus
	write.xlsx(selected_genes_with_degs, file = output_xlsx,sheetName = "DEGs_mm", append = FALSE)
	# genes DES m.musculus orthologs human
	write.xlsx(ortholog_table_res, file = output_xlsx,sheetName="DEGs_mm_ort_hs", append=TRUE)
	# genes obtained after feature selection
	if(directory_exp!="Comp8_C1A_vs_C1B"){
	write.xlsx(genes_to_save_FS, file = output_xlsx,sheetName="HS_FS", append=TRUE)
	}

	#
	# Prepare the data for tumour map analysis
	#

	mod2df<-data.frame(ID=names(mod2_classification),clusters_mclust=as.numeric(mod2_classification))
	mod2df$clusters_mclust<-as.character(paste("cluster",mod2df$clusters_mclust,sep="_"))
	
	ann_tot_for_tumap<-merge(ann3,mod2df,by.x="Samples",by.y="ID")[,-2]
	
	#
	#get only the samples IDH positive in GBM + LGG
	#	
	ann_idh_pos<-ann3[which(ann3$IDH1==1|ann3$IDH2==2),]
	ann_idh_pos2<-merge(ann_idh_pos,mod2df,by.x="Samples",by.y="ID")
	matrix_onlyIDH<-matrix_to_plot_hs_nocombat2[,which(colnames(matrix_to_plot_hs_nocombat2)%in%ann_idh_pos2[,1])]
	
	#
	# get only the samples IDH positive in LGG
	#	
	ann_idh1_lgg<-ann3[ann3$type=="LGG" & (ann3$IDH1 ==1 | ann3$IDH2 ==1 ),]
	ann_idh_pos2_lgg<-merge(ann_idh1_lgg,mod2df,by.x="Samples",by.y="ID")
	matrix_onlyIDH1_lgg<-matrix_to_plot_hs_nocombat2[,which(colnames(matrix_to_plot_hs_nocombat2)%in%ann_idh_pos2_lgg[,1])]
	
	#
	# Analysis methylation data
	#
  
	
	run_meth=TRUE

	# check if there are the data available from methylation
	if(run_meth==TRUE){

	genes_degs<-rownames(matrix_to_plot_hs_nocombat2)
	feature_selection_genes<-genes_to_save_FS
  
	tab_cpg_mouse<-read.delim(file="/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/regions_dmrs_adi_pval0.05_lfc1.txt")
	genes_cpg_mouse<-tab_cpg_mouse$annot.symbol
	
	ortholog_table_mgi_cpg<-ortholog_table[which(ortholog_table[,5]%in%genes_cpg_mouse),]
	#
	# The definitive version of biomaRt is 102 -- the default, to change the version use the parameter "version" December17 Ensembl Genes 102 ENSEMBL_MART_ENSEMBL
	#
	ortholog_table_ensembl_cpg<-unique(convertMouseGeneList2(genes_cpg_mouse))
	ortholog_table_res_meth<-integrated_conversion_mouse_human(mgi=ortholog_table_mgi_cpg,ensembl=ortholog_table_ensembl_cpg)
	
	CpG_island<-ortholog_table_res_meth[,1]
	  
	# sampleinfo_meth<-data.frame(ID=c("genes_degs","genes_degs","genes_degs","genes_degs","genes_degs","feature_selection_genes","CpG_island"),
	#            Rdata1=c(paste(directory_exp,pval,lfc,"meth.LGG_R132H.RData",sep="."),
	#                     paste(directory_exp,pval,lfc,"GBM_LGG_TRIPLE_MUT.RData",sep="."),
	#                     paste(directory_exp,pval,lfc,"meth.RData",sep="."),
	#                     paste(directory_exp,pval,lfc,"meth.LGG.RData",sep="."),
	#                     paste(directory_exp,pval,lfc,"meth.R132H.RData",sep="."),
	#                     paste(directory_exp,pval,lfc,"mclust",mkcluster,"meth.FS.RData",sep="."),
	#                     paste(directory_exp,pval,lfc,"mclust",mkcluster,"meth.FS.CpG.RData",sep=".")),
	#            output1=c(paste(directory_exp,pval,lfc,"meth.LGG_R132H.png",sep="."),
	#                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".GBM_LGG_TRIPLEMUT.png",sep=""),
	#                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".ge.meth.nocombat.png",sep=""),
	#                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".ge.meth.nocombat.LGG.png",sep=""),
	#                      paste(directory_exp,pval,lfc,"meth.R132H.png",sep="."),
	#                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".ge.meth.nocombat.FS.png",sep=""),
	#                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".ge.meth.nocombat.CpG.png",sep="")),
	#            output2=c(paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".meth.LGG_R132H.png",sep=""),
	#                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".only.METH.GBM_LGG_TRIPLEMUT.png",sep=""),
	#                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".ONLY.METH.nocombat.png",sep=""),
	#                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".ONLY.METH.nocombat.LGG.png",sep=""),
	#                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".meth.R132H.png",sep=""),
	#                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".ONLY.METH.nocombat.FS.png",sep=""),
	#                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".",mkcluster,".ONLY.METH.nocombat.CpG.png",sep="")),
	#            matrix_to_use=c("matrix_onlyR132H_lgg","matrix_TRIPLEMUT","matrix_to_plot_hs_nocombat2","matrix_to_plot_lgg","matrix_onlyR132H","matrix_to_plot_fs","matrix_to_plot_fs"),
	#            genes_to_use=c("genes_degs","genes_degs","genes_degs","genes_degs","genes_degs","feature_selection_genes","CpG_island"),
	#            qthr=c(20,20,20,20,20,20,20),stringsAsFactors = F)
	
	sampleinfo_meth<-data.frame(ID=c("genes_degs","genes_degs"),
	                            Rdata1=c(paste(directory_exp,pval,lfc,"meth.LGG_IDH_pos.Aug.RData",sep="."),
	                                     paste(directory_exp,pval,lfc,"meth.GBM_LGG_IDH_pos.Aug.RData",sep=".")),
	                            output1=c(paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),"ge.meth.LGG_IDH_pos.Aug.png",sep="."),
	                                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),"ge.meth.GBM_LGG_IDH_pos.Aug.png",sep=".")),
	                            output2=c(paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".meth.LGG_IDH_pos.Aug.png",sep=""),
	                                      paste(paste(directory_exp,"method.fc",pval,lfc,sep="."),".meth.GBM_LGG_IDH_pos.Aug.png",sep="")),
	                            matrix_to_use=c("matrix_onlyIDH1_lgg","matrix_onlyIDH"),
	                            genes_to_use=c("genes_degs","genes_degs"),
	                            qthr=c(20,20),stringsAsFactors = F)
	
	
	
	# December  15/12/2020: run the analysis only using the DEGs
	sampleinfo_meth<-sampleinfo_meth[sampleinfo_meth[,1]%in%c("genes_degs"),]
	           
	for(an_meth in 1:nrow(sampleinfo_meth)){
	
	print(an_meth)
	  
	type_experiment<-sampleinfo_meth[an_meth,1]  
	rdata_output<-sampleinfo_meth[an_meth,2]
	plot1_output<-as.character(sampleinfo_meth[an_meth,3])
	plot2_output<-as.character(sampleinfo_meth[an_meth,4])
	matrix_ge_to_use<-get(as.character(sampleinfo_meth[an_meth,5]))
	genes_to_use_meth<-get(as.character(sampleinfo_meth[an_meth,6]))
	qthr<-sampleinfo_meth[an_meth,7]
	
	ann_tot_for_tumap$IDH1_status<-rep(0,nrow(ann_tot_for_tumap))
	ann_tot_for_tumap[which(ann_tot_for_tumap$IDH1==1 | ann_tot_for_tumap$IDH2==1 & ann_tot_for_tumap$TP53==1 & ann_tot_for_tumap$ATRX==1),"IDH1_status"]<-1
		  
	library(data.table)
	library(doParallel)
	library(foreach)

	ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	tss <- getBM(attributes = c("transcription_start_site", "chromosome_name","transcript_start", "transcript_end", "strand",  "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"),filters = "external_gene_name", values = genes_to_use_meth,mart = ensembl)
	
	current_dir_gbm<-paste("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/GDCdata","TCGA-GBM","harmonized/DNA_Methylation/Methylation_Beta_Value",sep="/")
	meth_gbm<-create_matrix_meth(current_dir_gbm,genes_to_use_meth)
	#https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Methylation_LO_Pipeline/
  
	meth_gbm2<-meth_gbm[-which(meth_gbm[,7]%in%"."),]
	meth_gbm2<-meth_gbm2[-which(meth_gbm[,7]%in%"."),]
	meth_gbm3<-meth_gbm2[-which(is.na(meth_gbm2[,8])),]
	
	chromosome_cgi<-sapply(strsplit(gsub(meth_gbm3[,c(6)],30,pattern="CGI:",replacement=""),split=":"),"[[",1)
	start_cgi<-sapply(strsplit(sapply(strsplit(gsub(meth_gbm3[,c(6)],30,pattern="CGI:",replacement=""),split=":"),"[[",2),split="-"),"[[",1)
	end_cgi<-sapply(strsplit(sapply(strsplit(gsub(meth_gbm3[,c(6)],30,pattern="CGI:",replacement=""),split=":"),"[[",2),split="-"),"[[",2)
	cgi_df<-data.frame(chromosome_cgi,start_cgi=as.integer(start_cgi),end_cgi=as.integer(end_cgi),CGI=meth_gbm3[,c(6)],gene_symbol=meth_gbm3[,c(4)],stringsAsFactors = F)
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
  	meth_gbm4<-meth_gbm3[meth_gbm3[,6]%in%res_cgi_to_use,]

	print("create methylation matrix LGG")
	current_dir_lgg<-paste("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/GDCdata","TCGA-LGG","harmonized/DNA_Methylation/Methylation_Beta_Value",sep="/")
	
	meth_lgg<-create_matrix_meth(current_dir_lgg,genes_to_use_meth)
	meth_lgg2<-meth_lgg[-which(meth_lgg[,7]%in%"."),]
	meth_lgg2<-meth_lgg2[-which(meth_lgg[,7]%in%"."),]
	meth_lgg3<-meth_lgg2[-which(is.na(meth_lgg2[,8])),]
	
	chromosome_cgi<-sapply(strsplit(gsub(meth_lgg3[,c(6)],30,pattern="CGI:",replacement=""),split=":"),"[[",1)
	start_cgi<-sapply(strsplit(sapply(strsplit(gsub(meth_lgg3[,c(6)],30,pattern="CGI:",replacement=""),split=":"),"[[",2),split="-"),"[[",1)
	end_cgi<-sapply(strsplit(sapply(strsplit(gsub(meth_lgg3[,c(6)],30,pattern="CGI:",replacement=""),split=":"),"[[",2),split="-"),"[[",2)
	cgi_df<-data.frame(chromosome_cgi,start_cgi=as.integer(start_cgi),end_cgi=as.integer(end_cgi),CGI=meth_lgg3[,c(6)],gene_symbol=meth_lgg3[,c(4)],stringsAsFactors = F)
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
	meth_lgg4<-meth_lgg3[meth_lgg3[,6]%in%res_cgi_to_use,]
	
	#
	# Finally combine the methylation data of gbm and lgg
	#
	
	if(nrow(meth_lgg4)>nrow(meth_gbm4)){

	diff_to_correct<-nrow(meth_lgg4)-nrow(meth_gbm4)
	genes_gbm_meth_df<-data.frame(table(meth_gbm4[,4]))
	genes_lgg_meth_df<-data.frame(table(meth_lgg4[,4]))
	genes_gbm_lgg_all<-cbind(genes_gbm_meth_df,genes_lgg_meth_df)
	genes_to_correct<-genes_gbm_lgg_all[genes_gbm_lgg_all[,2]!=genes_gbm_lgg_all[,4],3]
	
	idx_genes_to_remove_randomly<-sample(which(meth_lgg4[,4]%in%genes_to_correct),1)
	meth_lgg4<-meth_lgg4[-idx_genes_to_remove_randomly,]

	}

	mat_meth<-cbind(meth_gbm4[,-c(c(1:3),c(5:7))],meth_lgg4[,-c(1:7)])
	colnames(mat_meth)<-gsub(unlist(lapply(strsplit(colnames(mat_meth),split="\\."),FUN=function(x){paste(x[1:3],collapse="-")})),pattern="\\.",replacement="-")
	        
 	setwd(output_dir)

	save(mat_meth,file=as.character(rdata_output))
	
	if(type_experiment=="genes_degs"){
	   
	plot_heatmap_hs_ge_meth(matrix_to_plot_hs=matrix_ge_to_use,
	 			mat_meth=mat_meth,
	 			ann_tot_for_tumap,
	 			output_hs=plot1_output,
	 			output_hs2=plot2_output,
	 			column_title="DEGs genes and associated meth islands",
	 			width=15,
	 			height=10,
	 			split=T,
	 			km=F,
	 			qthr=qthr,
	 			number_clusters=3,
	 			number_clusters_row=3,
	 			number_clusters_row_meth=3,
	 			how_split=mod2_classification)
	 
	 } else {
	   
	   plot_heatmap_hs_ge_meth(matrix_to_plot_hs=matrix_ge_to_use,
	                           mat_meth=mat_meth,
	                           ann_tot_for_tumap,
	                           output_hs=plot1_output,
	                           output_hs2=plot2_output,
	                           column_title="DEGs genes and associated meth islands",
	                           width=15,
	                           height=10,
	                           split=T,
	                           km=T,
	                           qthr=qthr,
	                           number_clusters=3,
	                           number_clusters_row=3,
	                           number_clusters_row_meth=3,
	                           how_split=mod2_classification)
	   
	   
	 }

	}
	
	}
	
	##########################################

	
	###############################################
	# Plot mouse results: all degs genes
 	###############################################	
  require(ggpubr)
	
	setwd(output_dir)
	
 	output_mm_boxplot<-paste("Mus_musculus",directory_exp,".method.fc",pval,lfc,".boxplot.sig.png",sep="")

	# png(output_mm_boxplot)
	# par(mar=c(15,3,1,1))
	# boxplot(log(matrix_to_plot_ms+1,2),las=2)
	# dev.off()
		
	output_mm<-paste("Mus_musculus",directory_exp,".method.fc",pval,lfc,".png",sep="")
	heatmapMouse(matrix_to_plot_ms,output_mm,ann_mm_clean)

  scores_mm<-CrossSpeciesComparison(matrix_to_plot_ms,ortholog_table_res,species="mm")
  ann_mm_clean_temp<-data.frame(ID=rownames(ann_mm_clean),ann_mm_clean)
  ann_mm_clean2<-merge(ann_mm_clean_temp,scores_mm,sort = F)
  
  output_mm<-paste("Mus_musculus",directory_exp,".method.fc",pval,lfc,".HISTscore.pdf",sep="")
  
  # pdf(output_mm)
  # 
  # phist<-gghistogram(ann_mm_clean2, x = "gsva_scores",
  #             add = "mean", rug = TRUE,
  #             color = "Condition", fill = "Condition",
  #             palette = c("#00AFBB", "#E7B800"),facet.by="Condition")
  # 
  # # print(phist)
  # 
  # pboxplot<-ggboxplot(ann_mm_clean2, "Condition", "gsva_scores",
  #           fill = "Condition", palette = c("#00AFBB", "#E7B800"))
  # 
  # print(pboxplot)
  # 
  # dev.off()
  
  # 
  ann_mm_clean_temp<-data.frame(ID=rownames(ann_mm_clean),ann_mm_clean)
  ann_mm_clean2<-merge(ann_mm_clean_temp,scores_mm,sort = F)
  # 
  output_mm<-paste("Mus_musculus",directory_exp,".method.fc",pval,lfc,".withSCORE.png",sep="")	        
  rownames(ann_mm_clean2)<-ann_mm_clean2[,1]
  heatmapMouse(matrix_to_plot_ms,output_mm,ann_mm_clean2[,-1])

	# if(feature_selection==1){
	# 
	# output_mm<-paste("Mus_musculus",directory_exp,".method.fc",pval,lfc,".png",sep="")
	# genes_to_save_FS_mm<-ortholog_table_res[ortholog_table_res[,1]%in%genes_to_save_FS,2]
	# matrix_to_plot_ms2<-matrix_to_plot_ms[rownames(matrix_to_plot_ms)%in%genes_to_save_FS_mm,]
	# colnames(matrix_to_plot_ms2)<-colnames(matrix_to_plot_ms)
	# 
	# output_txt<-paste("Mus_musculus",directory_exp,".method.fc",pval,lfc,".FS.",mkcluster,".txt",sep="")
	# write.table(genes_to_save_FS_mm,file=output_txt,sep="\t")
	# 
	# output_mm<-paste("Mus_musculus",directory_exp,".method.fc",pval,lfc,".FS.",mkcluster,".png",sep="")
	# heatmapMouse(matrix_to_plot_ms2,output_mm,ann_mm_clean)
	# 
	# }
	# 
	###############################################
	# Plot mouse results: ort genes with human
	############################################### 

	output_mm<-paste("Mus_musculus",directory_exp,".method.ort.fc",pval,lfc,".png",sep="")
	heatmapMouse(matrix_to_plot_ms_ort,output_mm,ann_mm_clean)

	setwd("/mnt/data/lab/gmt_data/data_brain/adi_february_2021")
	  
}	
