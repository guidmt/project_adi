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

prepareMatrixForCox<-function(mat_exp,annotation_file){
  
  matrix_for_gsva<-log(mat_exp+1,2) #convert in log scale for gsva
  
  matrix_for_mut_temp<-matrix_for_gsva[,which(colnames(matrix_for_gsva)%in%annotation_file[,1])] #get the ge samples that are in the annotation file
  nids<-data.frame(table(colnames(matrix_for_mut_temp))) #find the repeated samples
  samples_to_average<-nids[nids[,2]>1,1] 
  matrix_for_mut<-average_samples(matrix_for_mut_temp,samples_to_average) #collapse the expression values of the replicated samples
  
  return(matrix_for_mut)
  
}

# Define a function for gsva analysis
funGSVA<-function(file_gsva,X,min.sz=1,max.sz=2000,cut_off_gsva=0.5,cut_off_logfc=1,method="singscore",cox_genes=F,list_genes_cox=NULL){
  
  genes_for_gsva<-read.delim(file=file_gsva,stringsAsFactors=F)
  
  ovp_gsva_fs<-length(intersect(genes_for_gsva[,1],rownames(X)))
  
  print(paste("The overlap between C1AC1B and X is:",ovp_gsva_fs))
  
  ovp_gsva_fs2<-length(intersect(list_genes_cox,rownames(X)))
  
  print(paste("The overlap between COXPH genes and X is:",ovp_gsva_fs2))
  
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
system(paste("cp","/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/run_survival_analysis_final.R",xdate))

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

# sampleinfo<-read.delim(file="sampleinfo_C1vsC2.txt",header=F,stringsAsFactors=F,comment.char="#")

sampleinfo<-read.delim(file="sampleinfo_C1AvsC1B.txt",header=F,stringsAsFactors=F,comment.char="#")

RUN_CODE = TRUE

if(RUN_CODE == TRUE){
  
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
  # save_dir<-getwd()
  # system(paste("mkdir -p",paste(bd,directory_exp,sep="/")))
  
  setwd(paste("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/res_04_02_2021",paste(bd,directory_exp,sep="/"),sep="/"))
  output_dir<-getwd()
  
  #
  # Define the working directory
  #
  
  current_dir<-paste("/mnt/data/lab/gmt_data/data_brain",bd,directory_exp,sep="/")
  
  setwd(current_dir)
  
  input_file<-grep(dir(),pattern=".sigDE.csv",value=T)
  
  #
  # Here there is teh code used to select the genes in run_compare_general.R: just for double check
  #
  
  # degs_in_mouse<-read.csv(file=input_file,sep=",",check.names=F)
  # 
  # list_genes<-as.character(degs_in_mouse$Gene_name)
  # 
  # statistics_degs<-degs_in_mouse[,c(1:7)]
  # 
  # #	
  # # 4) Filter based on DEGs analysis 
  # #
  # 
  # #get the matrix with the gene-expression values
  # statistics_degs<-degs_in_mouse[,c(1:8)]
  # 
  # #select the differentially expressed genes
  # selected_genes_with_degs<-as.character(statistics_degs[statistics_degs$padj<pval & abs(statistics_degs$log2FoldChange)>=lfc,1])
  # 
  # #select the up and down regulated genes
  # up_genes<-statistics_degs[statistics_degs$padj<pval & statistics_degs$log2FoldChange>=lfc,1]
  # down_genes<-statistics_degs[statistics_degs$padj<pval & statistics_degs$log2FoldChange <= -c(lfc),1]
  # 
  # # save the DEGs in mm
  # output_degs<-paste(paste(directory_exp,pval,lfc,sep="_"),".txt",sep="")
  # 
  # output_txt<-paste(paste(output_dir,output_degs,sep="/"))
  # 
  # # write.table(selected_genes_with_degs,file=output_txt,sep="\t",row.names=T,quote=F)
  # 
  # #obtain the ortologs genes with hs
  # ortholog_table_mgi<-ortholog_table[which(ortholog_table[,5]%in%selected_genes_with_degs),]
  # #
  # # The definitive version of biomaRt is 102 -- the default, to change the version use the parameter "version"
  # #
  # print("WARNING:: The version of biomaRt used on 04/Feb/2021  (EU date) is 102, use convertMouseGeneList2 to use a old version ")
  # print("WARNING:: The version of biomaRt used on 04/Feb/2021  (EU date) is 102, use convertMouseGeneList2 to use a old version")
  # print("WARNING:: The version of biomaRt used on 04/Feb/2021  (EU date) is 102, use convertMouseGeneList2 to use a old version")
  # print("WARNING:: The version of biomaRt used on 04/Feb/2021  (EU date) is 102, use convertMouseGeneList2 to use a old version")
  # print("WARNING:: The version of biomaRt used on 04/Feb/2021  (EU date) is 102, use convertMouseGeneList2 to use a old version")
  # 
  # ortholog_table_ensembl<-unique(convertMouseGeneList2(selected_genes_with_degs))
  # 
  # ortholog_table_res<-integrated_conversion_mouse_human(mgi=ortholog_table_mgi,ensembl=ortholog_table_ensembl)
  # print("WARNING:: I am using the ensembl version 102, 04/Feb/2021. We locked the analysis using this repository.")
  # 
  output_degs2<-paste(paste(directory_exp,pval,lfc,sep="_"),".ort.txt",sep="")
  output_txt2<-paste(paste(output_dir,output_degs2,sep="/"))
  # 
  #
  # The file for GSVA analysis is the list of ort. genes in human - this analysis works only for C1A and C1B not C1 and C2
  #
  
  file_gsva<-output_txt2
  

  #
  # Create matrix h.sapiens --NOT-- processed with combat
  #
  
  col_to_select<-c(colnames(input_ge[,c(1:6)]),"Type")
  matrix_to_plot_hs_nocombat<-t(input_ge[,-c(which(colnames(input_ge)%in%col_to_select))])
  colnames(matrix_to_plot_hs_nocombat)<-unlist(lapply(X=strsplit(input_ge$Sample,split="\\-"),FUN=function(X)paste(X[1:3],collapse="-")))
  
  setwd(output_dir)
  
  system("rm all_results_coxph_adi.txt")

  #
  # Step 1: Load the clinical data of LGG and GBM with GDCquery_clinic
  #
  
      
  #get the samples with IDH1 + TP53 + ATRX mutations: N = 164 (There are 161 IDH1 and 3 IDH3)
  annotation_mut<-ann3[ann3$type%in%"LGG" & (ann3$IDH1==1 | ann3$IDH2==1)  & (ann3$TP53==1 & ann3$ATRX==1),]
  #get the samples with LGG and IDH1: N = 411
  annotation_only_idh<-ann3[ann3$type%in%"LGG" & (ann3$IDH1==1 | ann3$IDH2==1),]
  #get the samples that are LGG and GBM: N = 794
  annotation_gbm_lgg<-ann3
  #get the samples that are IDHmut-codel and IDHmut-non-codel:  414
  annotation_idh_codel_nocodel<-ann3[ann3$type%in%"LGG" & (ann3$IDH_codel_subtype=="IDHmut-codel" | ann3$IDH_codel_subtype=="IDHmut-non-codel"),]
  #get the samples that are LGG and IDHmut-codel:  169
  annotation_idh_codel<-ann3[ann3$type%in%"LGG" & (ann3$IDH_codel_subtype=="IDHmut-codel"),]
  #get the samples that are LGG and IDHmut-non-codel:  245
  annotation_idh_nocodel<-ann3[ann3$type%in%"LGG" & (ann3$IDH_codel_subtype=="IDHmut-non-codel"),]
  #get the samples that are GBM:  245
  annotation_gbm_only<-ann3[ann3$type%in%"GBM",]
      
  #
  # Step2: prepare the matrices with the gene-expression values for the analysis
  #
      
  matrices_for_gsva<-vector(mode="list",7)
  names(matrices_for_gsva)<-c("LGG_TRIPLE_MUT","only_LGG_IDH","GBM_LGG","IDH_CODEL_NOCODEL","IDH_CODEL","IDH_NOCODEL","only_GBM")
      
  list_annotations<-c("annotation_mut","annotation_only_idh","annotation_gbm_lgg","annotation_idh_codel_nocodel","annotation_idh_codel","annotation_idh_nocodel","annotation_gbm_only")
      
  matrix_for_gsva<-matrix_to_plot_hs_nocombat
      
  for(lafm in 1:length(list_annotations)){
        
      annotation_temp_for_cox<-get(list_annotations[lafm])  
      mat_temp_for_cox<-prepareMatrixForCox(matrix_for_gsva,annotation_temp_for_cox)
        
      matrices_for_gsva[[lafm]]<-mat_temp_for_cox
        
  }
  
  #
  # 2.1: prepare the clinical data with the overall survival
  #
  
  clin.lgg <- GDCquery_clinic("TCGA-LGG", "clinical")
  clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
  
  features_to_select<-c("submitter_id","days_to_death","days_to_last_follow_up","vital_status")
  
  clin.lgg2<-clin.lgg[,colnames(clin.lgg)%in%features_to_select]
  clin.gbm2<-clin.gbm[,colnames(clin.gbm)%in%features_to_select]
  
  #
  # Step3: Select the genes only using the samples that are IDH1/2 (MUT) + TP53 (MUT) + ATRX (MUT)
  #
  
  group_of_start<-names(matrices_for_gsva)[c(1,6)]
  
  ann_of_start<-c("annotation_mut","annotation_only_idh","annotation_gbm_lgg","annotation_idh_codel_nocodel","annotation_idh_codel","annotation_idh_nocodel","annotation_gbm_only")
  names(ann_of_start)<-names(matrices_for_gsva)
  
  res_list_cox<-vector(mode="list",length(group_of_start))
  names(res_list_cox)<-group_of_start
    
  for(an_gs in 1:length(group_of_start)){
    
    an_gs_string<-group_of_start[an_gs]
    # remove the folder of results if exist
    remove_folder_res<-paste(paste("rm -r",an_gs_string))
    
    system(remove_folder_res)
    
    create_folder_res<-paste(paste("mkdir",an_gs_string))
    system(create_folder_res)
    
    features_to_select<-c("submitter_id","days_to_death","days_to_last_follow_up","vital_status")
    #this clinical data contains 515 patients, but they are merged with the selected samples in the gene-expression (function coxphAnalysis)
    # TCGA_surv_coxph<-clin.lgg2[,colnames(clin.lgg2)%in%features_to_select]
    TCGA_surv_coxph_all<-rbind(clin.lgg2,clin.gbm2)
    
    # get the annotation data derived previosuly to analyze
    temp_ann_coxph<-get(ann_of_start[an_gs_string])
    
    #use only the samples in the annotation file 
    TCGA_surv_coxph<-TCGA_surv_coxph_all[which(TCGA_surv_coxph_all[,1]%in%temp_ann_coxph[,1]),]
    
    # the maximum cut-off of significance is 0.05
    resCox<-coxphAnalysis(TCGA_surv_coxph,matrices_for_gsva[[an_gs]],file_gsva,thr_cox=0.05)

    # mtx_for_gsva_coxgenes<-matrices_for_gsva[[an_gs]]
    # mtx_for_gsva_coxgenes2<-mtx_for_gsva_coxgenes[rownames(mtx_for_gsva_coxgenes)%in%resCox_genes,]
    
    # samplesc1ac1b<-gsub(sapply(strsplit(file_gsva,split="/"),"[[",13),pattern=".ort.txt",replacement="")
    
    res_list_cox[[an_gs]]<-resCox
      
  }
  
    
  #
  # Step4: - Use the genes of IDH1/2 (MUT) + TP53 (MUT) + ATRX (MUT) to see the results of the O.S. in other groups and different thr of p.values
  #        - Of course the procedure run the analysis also considering the group of the TRIPLE mutants
  
  for(get_cox in 1:length(res_list_cox)){
    
  res_cox_to_use<-res_list_cox[[get_cox]]
  names_cox_to_use<-names(res_list_cox[get_cox])
  
  create_folder_res_full_path<-paste(getwd(),names_cox_to_use,sep="/")
  
  system(paste("rm",paste(output_dir,"/",paste(names_cox_to_use,"all_results_coxph_adi.txt",sep="_"),sep="")))
  
  # for each group to analyse
  for(an_gs2 in 1:length(matrices_for_gsva)){
      
  # Check if the DEGs C1A and C1B are able to distinguish several biological groups (independently by coxph)
  scores_gsva_degs_genes<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs2]],method="gsva",max.sz=2000,cut_off_gsva=0)
  current_score_gsva_degs_genes<-setDT(scores_gsva_degs_genes[[2]])
  all_clinical_for_degs_genes<-rbind(clin.lgg2,clin.gbm2)
  
  input_os_degs_genes<-merge(all_clinical_for_degs_genes,current_score_gsva_degs_genes,by.x="submitter_id",by.y="ID")
  input_os_degs_genes$status_C1AC1B_new<-as.character(input_os_degs_genes$status_C1AC1B)
  
  input_os_degs_genes[which(input_os_degs_genes$global_score>=-0.5 & input_os_degs_genes$global_score<=0.5),"status_C1AC1B_new"]<-"uncategorized"
  input_os_degs_genes_rm<-input_os_degs_genes[-which(input_os_degs_genes$status_C1AC1B_new%in%"uncategorized"),]
  
  output_string<-paste(scores_gsva_degs_genes[[1]],sep="_")
  
  # Overall Survival create using all the genes DEGs between C1A and C1B and a cut-off ssgsea of 0
  p1<-plotSurvival(TCGA_surv=input_os_degs_genes,status="status_C1AC1B",title=paste(output_string,"custom thr. DEGs C1A/C1B"))
  # Overall Survival create using all the genes DEGs between C1A and C1B and a cut-off ssgsea of 0, with "uncategorized" samples
  p2<-plotSurvival(TCGA_surv=input_os_degs_genes,status="status_C1AC1B_new",title=paste("with uncategorized","custom thr. DEGs C1A/C1B"))
  # Overall Survival create using all the genes DEGs between C1A and C1B and a cut-off ssgsea of 0, withOUT "uncategorized" samples
  p3<-plotSurvival(TCGA_surv=input_os_degs_genes_rm,status="status_C1AC1B_new",title=paste("without uncategorized","custom thr. DEGs C1A/C1B"))
  
  # Overall Survival create using all the genes DEGs between C1A and C1B using automatic threshold
  p4<-plotSurvival2(TCGA_surv=input_os_degs_genes,score_to_use="global_score",suffix_output="aut.thr with UN")[[3]]
  p5<-plotSurvival2(TCGA_surv=input_os_degs_genes_rm,score_to_use="global_score",suffix_output="aut.thr without UN")[[3]]
  
  pdf(paste("ALL_DEGS_C1A_C1B",pval,lfc,names(matrices_for_gsva)[an_gs2],"pdf",sep="."),width=10,height = 18)
  
  survplotlist<-list()
  survplotlist[[1]]<-p1
  survplotlist[[2]]<-p2
  survplotlist[[3]]<-p3
  survplotlist[[4]]<-p4
  survplotlist[[5]]<-p5
  
  arrange_ggsurvplots(survplotlist,ncol=2,nrow=3,print=T,title=names(matrices_for_gsva)[[an_gs2]])
  
  dev.off()
  
  
      #define threhsolds  for the logRank analysis
      list_thr_cox<-c(0.01,0.001,0.05)
      
      for(thr_cox in list_thr_cox){
        
        # Use the genes selected in Step3
        resCox_genes<-data.frame(res_cox_to_use[res_cox_to_use$LogRank<thr_cox,"Variable"])[,1]
        resCox_genes_select<-res_cox_to_use[which(res_cox_to_use$LogRank<thr_cox),]
        
        print(paste("i am using the genes selected in step3:",length(resCox_genes)))
        
        # mtx_for_gsva_coxgenes<-matrices_for_gsva[[an_gs2]]
        # mtx_for_gsva_coxgenes2<-mtx_for_gsva_coxgenes[rownames(mtx_for_gsva_coxgenes)%in%resCox_genes,]
        
        samplesc1ac1b<-gsub(sapply(strsplit(file_gsva,split="/"),"[[",13),pattern=".ort.txt",replacement="")
        
        # define the type of gene list to use
        list_degs_cox<-c("up_dw","up","dw")
        
        for(type_degs_cox in list_degs_cox){
          
          print(paste("Analysis:",names(matrices_for_gsva)[[an_gs2]],"-",thr_cox,"-",type_degs_cox))
          
          if(type_degs_cox == "up_dw"){
            
            print(type_degs_cox)
            
            scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs2]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes,cut_off_gsva=0)
            
          } 
          
          if(type_degs_cox == "up"){
            
            print(type_degs_cox)
            
            # select only up-regulated genes 
            resCox_genes_up<-as.character(resCox_genes_select[resCox_genes_select$status %in% "up",]$Variable)
            
            #if there are not sufficient up regulated genes from the cox analysis use directly the up regulated genes from the DEGs analysis
            if(length(resCox_genes_up)==0){
              
              resCox_genes_up<-as.character(ortholog_table_res[ortholog_table_res[,3]%in%"up",1])
              
              scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs2]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_up,cut_off_gsva=0)
              
            } else{
              
              # use the up regulated genes from coxph   
              resCox_genes_up<-as.character(resCox_genes_select[resCox_genes_select$status %in% "up",]$Variable)
              
              scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs2]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_up,cut_off_gsva=0)
              
            }
            
          }
          
          if(type_degs_cox == "dw"){
            
            print(type_degs_cox)
            
            # select only down-regulated genes
            resCox_genes_dw<-as.character(resCox_genes_select[resCox_genes_select$status %in% "down",]$Variable)
            
            #if there are not sufficient down regulated genes from the cox analysis use directly the down regulated genes from the DEGs analysis
            if(length(resCox_genes_dw)==0){
              
              resCox_genes_dw<-as.character(ortholog_table_res[ortholog_table_res[,3]%in%"down",1])
              
              scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs2]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_dw,cut_off_gsva=0)
              
            } else {
              
              # use the down regulated genes from coxph   
              resCox_genes_dw<-as.character(resCox_genes_select[resCox_genes_select$status %in% "down",]$Variable)
              
              scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs2]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_dw,cut_off_gsva=0)
              
            } 
            
          } 
          
          current_score_gsva<-scores_gsva[[2]]
          
          output_string<-paste(scores_gsva[[1]],names(matrices_for_gsva)[[an_gs2]],sep="_")
          
          output_pdf<- paste(create_folder_res_full_path,"/",paste(paste("OS",sep="."),samplesc1ac1b,output_string,"coxph_thrcox",thr_cox,type_degs_cox,"pdf",sep="."),sep="")
          
          
          output_cox<-paste(create_folder_res_full_path,"/",paste(paste("COX_RESULTS_thr_cox",thr_cox,sep="."),samplesc1ac1b,output_string,"coxph.txt",sep="."),sep="")
          write.table(resCox_genes_select,file=output_cox,sep="\t",row.names=F,quote=F)
          # 
          # output_os_txt<- paste(paste("OS",sep="."),samplesc1ac1b,output_string,"coxph_thrcox",thr_cox,type_degs_cox,"txt",sep=".")
          # write.table(current_score_gsva,file=output_os_txt,sep="\t",row.names=F,quote=F)
          # 
          # outpdf2<-paste(getwd(),output_pdf,sep="/")
          # outtxt2<-paste(getwd(),output_cox,sep="/")
          
          outpdf2<-output_pdf
          outtxt2<-output_cox
          
          line=paste(names(matrices_for_gsva)[[an_gs2]],outpdf2,outtxt2,thr_cox,type_degs_cox,sep="\t")
          
          write(line,file=paste(names_cox_to_use,"all_results_coxph_adi.txt",sep="_"),append=TRUE)
          

          TCGA_surv_coxph2<-TCGA_surv_coxph_all[which(TCGA_surv_coxph_all[,1]%in%current_score_gsva[,1]),]
          
          input_os_from_cox<-merge(TCGA_surv_coxph2,current_score_gsva,by.x="submitter_id",by.y="ID",all.x=F,all.y=T)
          
          #
          # There are samples in which the GSVA scores is between -0.5/0.5, we define this as a third group
          #
          
          input_os_from_cox_for_bordeline<-input_os_from_cox
          input_os_from_cox_for_bordeline$status_C1AC1B_new<-as.character(input_os_from_cox_for_bordeline$status_C1AC1B)
          
          input_os_from_cox_for_bordeline[which(input_os_from_cox_for_bordeline$global_score>=-0.5 & input_os_from_cox_for_bordeline$global_score<=0.5),"status_C1AC1B_new"]<-"uncategorized"
          input_os_from_cox_for_bordeline_rm<-input_os_from_cox_for_bordeline[-which(input_os_from_cox_for_bordeline$status_C1AC1B_new%in%"uncategorized"),]
          
          # 1: make the survival with the three group
          p1<-plotSurvival(TCGA=input_os_from_cox_for_bordeline,status="status_C1AC1B_new",title=paste(output_string,"custom thr",type_degs_cox))
          # 2: make the survival without the uncategorized group
          
          if(nrow(input_os_from_cox_for_bordeline_rm)>10){
            
          p2<-plotSurvival(TCGA=input_os_from_cox_for_bordeline_rm,status="status_C1AC1B_new",title=paste(output_string,"custom thr",type_degs_cox))
          # 3: get the chart of the automatic split of the data without the uncategorized group, and automatic cut-off  
          p3<-plotSurvival2(TCGA_surv=input_os_from_cox_for_bordeline_rm,score_to_use="global_score",suffix_output=paste(type_degs_cox))[[2]]
          # 4: get the chart of the survival of the data without the uncategorized group, and automatic cut-off  
          p4<-plotSurvival2(TCGA_surv=input_os_from_cox_for_bordeline_rm,score_to_use="global_score",suffix_output=paste(type_degs_cox))[[3]]
          
          }
          
          outpdf3<- paste(create_folder_res_full_path,"/",paste(paste("OS",sep="."),samplesc1ac1b,names(matrices_for_gsva)[[an_gs2]],"with_UNCATEGORIZED",thr_cox,type_degs_cox,"pdf",sep="."),sep="")
          
          pdf(outpdf3,width=15)
          survplotlist<-list()
          survplotlist[[1]]<-p1
          if(nrow(input_os_from_cox_for_bordeline_rm)>1){
          survplotlist[[2]]<-p2
          survplotlist[[3]]<-p4
          }
          arrange_ggsurvplots(survplotlist,ncol=3,nrow=1,print=T)
          dev.off()
          
          pdf(output_pdf)
          
          #
          #  Use the function plotSurvival and column status_C1AC1B to plot the KM
          #
          plotSurvival(TCGA_surv=input_os_from_cox,status="status_C1AC1B",title=paste(output_string,"custom thr",type_degs_cox))
          
          #
          # Use the function plotSurvival2 and column global score to determine the optimal cut-off point and plot the KM
          #
          
          # res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="global_score",suffix_output=paste(output_string,"up_and_down"))
          res_os<-plotSurvival2(TCGA_surv=input_os_from_cox,score_to_use="global_score",suffix_output=paste("gsva_aut.thr.",type_degs_cox))[[1]]
          
          res_os2<-unique(res_os[,c(1,4)])
          colnames(res_os2)[c(1:2)]<-c("Samples","automatic_classification")
          
          res_os3<-merge(res_os2,unique(input_os_from_cox[,c(1,5,6,7,8)]),by.x="Samples",by.y="submitter_id")
          
          output_os_txt<- paste(create_folder_res_full_path,"/",paste(paste("OS",sep="."),samplesc1ac1b,output_string,"coxph_thrcox",thr_cox,type_degs_cox,"txt",sep="."),sep="")
          write.table(res_os3,file=output_os_txt,sep="\t",row.names=F,quote=F)
          
          dev.off()
          
        } # end loop for up/dw/up_and_down
        
      } # end loop for logRank thresholds
      
  } # end loop for groups
  
  } # end loop for results cox
  
  
  #
  # Create reports survival analysis
  #
  
  #Plot results TRIPLE MUT
  output_file<-paste("TRIPLE_MUT_report_survival_analysis_coxph_logfc=",lfc,"qvalue=",pval,"html",sep=".")
  rmarkdown::render('/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/res_04_02_2021/template_OS_Adi.TRIPLEMUT.Rmd', output_file = output_file)
  
  #Plot results IDH-NOCODEL
  output_file<-paste("IDH_NOCODEL_report_survival_analysis_coxph_logfc=",lfc,"qvalue=",pval,"html",sep=".")
  rmarkdown::render('/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/res_04_02_2021/template_OS_Adi.IDH_NOCODEL.Rmd', output_file = output_file)
  
  setwd("/mnt/data/lab/gmt_data/data_brain/adi_february_2021")
  
}	# end loop for sampleinfo

}
