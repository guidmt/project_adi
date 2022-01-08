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
    # Step3: download the clinical data
    #
    
    clin.lgg <- GDCquery_clinic("TCGA-LGG", "clinical")
    clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
    
    features_to_select<-c("submitter_id","days_to_death","days_to_last_follow_up","vital_status")
    
    clin.lgg2<-clin.lgg[,colnames(clin.lgg)%in%features_to_select]
    clin.gbm2<-clin.gbm[,colnames(clin.gbm)%in%features_to_select]
    
    TCGA_surv_all<-rbind(clin.lgg2,clin.gbm2)
    
    #
    # Step3: Scoring using Buffa based approach
    #

    degs_genes<-read.delim(file=file_gsva)
    
    up_c1a<-degs_genes[which(degs_genes$status%in%"up"),]
    up_c1b<-degs_genes[which(degs_genes$status%in%"down"),]
  
    gene_sets_list<-vector(mode="list",2)
    
    gene_sets_list[[1]]<-up_c1a
    gene_sets_list[[2]]<-up_c1b
    
    for(ib in 1:length(matrices_for_gsva)){
    
    current_matrix_ib<-matrices_for_gsva[[ib]]
    current_group<-names(matrices_for_gsva)[ib]
    
    results_scoring<-vector(mode="list",2)
    
    #for each signature
    for(scorb in 1:length(gene_sets_list)){
      
      # select the gene-expression data of the signature
      gene_set_hypoxia<-gene_sets_list[[scorb]][,1]
      temp_matrix_gs<-current_matrix_ib[which(rownames(current_matrix_ib)%in%gene_set_hypoxia),]
  
      # create a temporary matrix
      scores_by_gs<-vector(mode="list",nrow(temp_matrix_gs))
      
      # for each gene in the hypoxia signature
      for(gene_to_score in 1:nrow(temp_matrix_gs)){
        
        # select expression values for a gene
        temp_matrix_gs2<-as.numeric(temp_matrix_gs[gene_to_score,])
        
        # compute quantile 
        quantile_gene<-quantile(temp_matrix_gs2)
        
        # discretize values scores https://www.thelancet.com/journals/lanonc/article/PIIS1470-2045(14)71021-6/fulltext
        res_buffa_scoring<-ifelse(temp_matrix_gs2>quantile_gene[3],1,-1)
        
        scores_by_gs[[gene_to_score]]<-res_buffa_scoring
        
      }
      
      scores_current_gs<-apply(do.call(rbind,scores_by_gs),2,sum)
      names(scores_current_gs)<-colnames(temp_matrix_gs)
        
      results_scoring[[scorb]]<-scores_current_gs
      
      }
      
    df_score_c1a<-data.frame(ID=names(results_scoring[[1]]),scoreC1A=results_scoring[[1]])
    df_score_c1b<-data.frame(ID=names(results_scoring[[2]]),scoreC1B=results_scoring[[2]])
    
    df_tot<-cbind(df_score_c1a,df_score_c1b)
    df_tot$status_final<-rep("none",nrow(df_tot))
    
    df_tot[which(df_tot$scoreC1A<0 & df_tot$scoreC1B <0),"status_final"]<-"down_in_both_groups"
    df_tot[which(df_tot$scoreC1A>0 & df_tot$scoreC1B <0),"status_final"]<-"up_C1A"
    df_tot[which((df_tot$scoreC1A>df_tot$scoreC1B) & (df_tot$scoreC1A>0 & df_tot$scoreC1B>0)),"status_final"]<-"up_C1A"
    df_tot[which(df_tot$scoreC1A<0 & df_tot$scoreC1B >0),"status_final"]<-"up_C1B"
    df_tot[which((df_tot$scoreC1B>df_tot$scoreC1A) & (df_tot$scoreC1B>0 & df_tot$scoreC1A>0)),"status_final"]<-"up_C1B"
    df_tot[,1]<-as.character(df_tot[,1])
    colnames(df_tot)[1]<-"id_ref"
    
    TCGA_surv_select_groups<-TCGA_surv_all[which(TCGA_surv_all[,1]%in%df_tot[,1]),]
    TCGA_surv_select_groups2<-merge(TCGA_surv_select_groups,df_tot,by.x="submitter_id",by.y="id_ref")
    
    pdf(paste(current_group,"_Buffa_scoring.pdf",sep=""),width=12)
    
    p1<-plotSurvival(TCGA_surv=TCGA_surv_select_groups2,status="status_final",title=paste(current_group,"Buffa scoring"))
    
    p2<-plotSurvival(TCGA_surv=TCGA_surv_select_groups2[-which(TCGA_surv_select_groups2$status_final%in%"down_in_both_groups"),],status="status_final",title=paste(current_group,"Buffa scoring"))

    survplotlist<-list()
    survplotlist[[1]]<-p1
    survplotlist[[2]]<-p2

    arrange_ggsurvplots(survplotlist,ncol=2,nrow=1,print=T,title=current_group)

    dev.off()

    }
    
  }