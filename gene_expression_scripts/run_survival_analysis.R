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
system(paste("cp","/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/run_survival_analysis.R",xdate))

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
  
  annotation_all<-ann3[,c("Samples","type","IDH1","IDH1_R132H","IDH1_other","TP53","ATRX","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]
  
    if(do_gsva==1){
      
      #get thhe samples with IDH1 + TP53 + ATRX mutations
      annotation_mut<-annotation_all[annotation_all$type%in%"LGG" & (annotation_all$IDH1==1 & annotation_all$TP53==1 & annotation_all$ATRX==1),]
      #get thhe samples with LGG and IDH1
      annotation_only_idh<-annotation_all[annotation_all$type%in%"LGG" & (annotation_all$IDH1),]
      annotation_gbm_lgg<-annotation_all
      annotation_idh_codel_nocodel<-annotation_all[annotation_all$type%in%"LGG" & (annotation_all$IDH_codel_subtype=="IDHmut-codel" | annotation_all$IDH_codel_subtype=="IDHmut-non-codel"),]
      annotation_idh_codel<-annotation_all[annotation_all$type%in%"LGG" & (annotation_all$IDH_codel_subtype=="IDHmut-codel"),]
      annotation_idh_nocodel<-annotation_all[annotation_all$type%in%"LGG" & (annotation_all$IDH_codel_subtype=="IDHmut-non-codel"),]
      annotation_gbm_only<-annotation_all[annotation_all$type%in%"GBM",]
      
      #
      # Define a function to prepare the input for coxph
      #
      
      prepareMatrixForCox<-function(mat_exp,annotation_file){
        
        matrix_for_gsva<-log(mat_exp+1,2) #convert in log scale for gsva
        
        matrix_for_mut_temp<-matrix_for_gsva[,which(colnames(matrix_for_gsva)%in%annotation_file[,1])] #get the ge samples that are in the annotation file
        nids<-data.frame(table(colnames(matrix_for_mut_temp))) #find the repeated samples
        samples_to_average<-nids[nids[,2]>1,1] 
        matrix_for_mut<-average_samples(matrix_for_mut_temp,samples_to_average) #collapse the expression values of the replicated samples
        
        return(matrix_for_mut)
        
      }
      
      matrices_for_gsva<-vector(mode="list",7)
      names(matrices_for_gsva)<-c("LGG_TRIPLE_MUT","only_LGG_IDH","GBM_LGG","IDH_CODEL_NOCODEL","IDH_CODEL","IDH_NOCODEL","only_GBM")
      
      list_annotations<-c("annotation_mut","annotation_only_idh","annotation_gbm_lgg","annotation_idh_codel_nocodel","annotation_idh_codel","annotation_idh_nocodel","annotation_gbm_only")
      
      matrix_for_gsva<-matrix_to_plot_hs_nocombat
      
      for(lafm in 1:length(list_annotations)){
        
        annotation_temp_for_cox<-get(list_annotations[lafm])  
        mat_temp_for_cox<-prepareMatrixForCox(matrix_for_gsva,annotation_temp_for_cox)
        
        matrices_for_gsva[[lafm]]<-mat_temp_for_cox
        
      }
      

      #Download the clinical data for LGG and GBM
      clin.lgg <- GDCquery_clinic("TCGA-LGG", "clinical")
      clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
      
      #Download the important columns for the analysis
      features_to_select<-c("submitter_id","days_to_death","days_to_last_follow_up","vital_status")
      clin.lgg2<-clin.lgg[,colnames(clin.lgg)%in%features_to_select]
      clin.gbm2<-clin.gbm[,colnames(clin.gbm)%in%features_to_select]
      
      for(an_gs in 1:length(matrices_for_gsva)){
        
        # read the file with the differentially expressed genes           
        samplesc1ac1b<-gsub(sapply(strsplit(file_gsva,split="/"),"[[",13),pattern=".ort.txt",replacement="")
        
        # use the DEGs to compute a score with gsva
        scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",max.sz=2000,cut_off_gsva=0)
        
        # extract the name of the experiment
        current_score_gsva<-setDT(scores_gsva[[2]])
        
        output_string<-paste(scores_gsva[[1]],names(matrices_for_gsva)[[an_gs]],sep="_")
        
        output_pdf<- paste(paste("OS",sep="."),samplesc1ac1b,output_string,"pdf",sep=".")
        
        output_os_txt<- paste(paste("OS",sep="."),samplesc1ac1b,output_string,"txt",sep=".")
        write.table(scores_gsva,file=output_os_txt,sep="\t",row.names=F,quote=F)
        
        pdf(output_pdf)
        
        # put toghether the clinical data
        input_os_gsva_all<-rbind(clin.lgg2,clin.gbm2)
        
        # input_os_gsva_lgg<-merge(clin.lgg2,current_score_gsva,by.x="submitter_id",by.y="ID")
        
        # use only the clinical data selected available in the matrix with the scores
        input_os_gsva_curr_group<-merge(input_os_gsva_all,current_score_gsva,by.x="submitter_id",by.y="ID")
        
        #
        # os with cut-off approach 
        #
        
        # plotSurvival(TCGA_surv=input_os_gsva_lgg,status="status_C1AC1B",title=paste(output_string,"custom thr. up_and_down"))
        plotSurvival(TCGA_surv=input_os_gsva_curr_group,status="status_C1AC1B",title=paste(output_string,"custom thr. up_and_down"))
        
        #
        # os with best cut-point
        #
        
        
        # res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="global_score",suffix_output=paste(output_string,"up_and_down"))
        res_os<-plotSurvival2(TCGA_surv=input_os_gsva_curr_group,score_to_use="global_score",suffix_output=paste(output_string,"up_and_down"))
        
        #res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="down_score",suffix_output=paste(output_string,"only_down"))
        
        #res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="up_score",suffix_output=paste(output_string,"only_up"))
        
        dev.off()
      }
      
    }
    
    #
    # Survival analysis coxph 
    #
    system("rm all_results_coxph_adi.txt")
    
    clin.lgg <- GDCquery_clinic("TCGA-LGG", "clinical")
    clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
    
    clin.lgg2<-clin.lgg[,colnames(clin.lgg)%in%features_to_select]
    clin.gbm2<-clin.gbm[,colnames(clin.gbm)%in%features_to_select]
    
    for(an_gs in 1:length(matrices_for_gsva[[1]])){
      
      features_to_select<-c("submitter_id","days_to_death","days_to_last_follow_up","vital_status")
      #this clinical data contains 515 patients, but they are merged with the selected samples in the gene-expression (function coxphAnalysis)
      # TCGA_surv_coxph<-clin.lgg2[,colnames(clin.lgg2)%in%features_to_select]
      TCGA_surv_coxph_all<-rbind(clin.lgg2,clin.gbm2)
      
      # get the annotation data derived previosuly to analyze
      temp_ann_coxph<-get(list_annotations[an_gs])
      
      #use only the samples in the annotation file 
      TCGA_surv_coxph<-TCGA_surv_coxph_all[which(TCGA_surv_coxph_all[,1]%in%temp_ann_coxph[,1]),]
      
      #define threhsolds  for the logRank analysis
      list_thr_cox<-c(0.01,0.001,0.05)
      
      for(thr_cox in list_thr_cox){
        
        # find the important genes with Coxph
        resCox<-coxphAnalysis(TCGA_surv_coxph,matrices_for_gsva[[an_gs]],file_gsva,thr_cox=thr_cox)
        resCox_genes<-resCox$Variable
        
        mtx_for_gsva_coxgenes<-matrices_for_gsva[[an_gs]]
        mtx_for_gsva_coxgenes2<-mtx_for_gsva_coxgenes[rownames(mtx_for_gsva_coxgenes)%in%resCox_genes,]
        
        samplesc1ac1b<-gsub(sapply(strsplit(file_gsva,split="/"),"[[",13),pattern=".ort.txt",replacement="")
        
        # define the type of gene list to use
        list_degs_cox<-c("up_dw","up","dw")
        
        all_res_cox2<-data.frame()
        
        for(type_degs_cox in list_degs_cox){
          
          print(paste("Analysis:",names(matrices_for_gsva)[[an_gs]],"-",thr_cox,"-",type_degs_cox))
          
          if(type_degs_cox == "up_dw"){
            
            print(type_degs_cox)
            
            scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes,cut_off_gsva=0)
            
          } 
          
          if(type_degs_cox == "up"){
            
            print(type_degs_cox)
            
            # select only up-regulated genes 
            resCox_genes_up<-as.character(resCox[resCox$status %in% "up",]$Variable)
            
            #if there are not sufficient up regulated genes from the cox analysis use directly the up regulated genes from the DEGs analysis
            if(length(resCox_genes_up)==0){
              
              resCox_genes_up<-as.character(ortholog_table_res[ortholog_table_res[,3]%in%"up",1])
              
              scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_up,cut_off_gsva=0)
              
            } else{
              
              # use the up regulated genes from coxph   
              resCox_genes_up<-as.character(resCox[resCox$status %in% "up",]$Variable)
              
              scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_up,cut_off_gsva=0)
              
            }
            
          }
          
          if(type_degs_cox == "dw"){
            
            print(type_degs_cox)
            
            # select only down-regulated genes
            resCox_genes_dw<-as.character(resCox[resCox$status %in% "dw",]$Variable)
            
            #if there are not sufficient down regulated genes from the cox analysis use directly the down regulated genes from the DEGs analysis
            if(length(resCox_genes_dw)==0){
              
              resCox_genes_dw<-as.character(ortholog_table_res[ortholog_table_res[,3]%in%"down",1])
              
              scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_dw,cut_off_gsva=0)
              
            } else {
              
              # use the down regulated genes from coxph   
              resCox_genes_dw<-as.character(resCox[resCox$status %in% "down",]$Variable)
              
              scores_gsva<-funGSVA(file_gsva=file_gsva,X=matrices_for_gsva[[an_gs]],method="gsva",cox_genes=T,list_genes_cox=resCox_genes_dw,cut_off_gsva=0)
              
            } 
            
          } 
          
          current_score_gsva<-scores_gsva[[2]]
          
          output_string<-paste(scores_gsva[[1]],names(matrices_for_gsva)[[an_gs]],sep="_")
          
          output_pdf<- paste(paste("OS",sep="."),samplesc1ac1b,output_string,"coxph_thrcox",thr_cox,type_degs_cox,"pdf",sep=".")
          
          
          output_cox<-paste(paste("COX_RESULTS_thr_cox",thr_cox,sep="."),samplesc1ac1b,output_string,"coxph.txt",sep=".")
          write.table(resCox,file=output_cox,sep="\t",row.names=F,quote=F)
          
          output_os_txt<- paste(paste("OS",sep="."),samplesc1ac1b,output_string,"coxph_thrcox",thr_cox,type_degs_cox,"txt",sep=".")
          write.table(current_score_gsva,file=output_os_txt,sep="\t",row.names=F,quote=F)
          
          outpdf2<-paste(getwd(),output_pdf,sep="/")
          outtxt2<-paste(getwd(),output_cox,sep="/")
          
          line=paste(names(matrices_for_gsva)[[an_gs]],outpdf2,outtxt2,thr_cox,type_degs_cox,sep="\t")
          
          write(line,file="all_results_coxph_adi.txt",append=TRUE)
          
          pdf(output_pdf)
          
          #
          # Prepare clinical data
          #
          
          # input_os_gsva_lgg<-merge(clin.lgg2,current_score_gsva,by.x="submitter_id",by.y="ID")
          input_os_from_cox<-merge(TCGA_surv_coxph,current_score_gsva,by.x="submitter_id",by.y="ID",all.x=F,all.y=T)
          
          #
          # os with cut-off approach - up and down genes 
          #
          
          # plotSurvival(TCGA_surv=input_os_gsva_lgg,status="status_C1AC1B",title=paste(output_string,"custom thr. up_and_down"))
          plotSurvival(TCGA_surv=input_os_from_cox,status="status_C1AC1B",title=paste(output_string,"custom thr",type_degs_cox))
          
          #
          # os with best cut-point - up and down
          #
          
          # res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="global_score",suffix_output=paste(output_string,"up_and_down"))
          res_os<-plotSurvival2(TCGA_surv=input_os_from_cox,score_to_use="global_score",suffix_output=paste(output_string,type_degs_cox))
          
          #res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="down_score",suffix_output=paste(output_string,"only_down"))
          
          #res_os<-plotSurvival2(TCGA_surv=input_os_gsva_lgg,score_to_use="up_score",suffix_output=paste(output_string,"only_up"))
          
          dev.off()
          
        }
        
      }
      
    }
  
  #
  # Create reports survival analysis
  #
  
  output_file<-paste("report_survival_analysis_coxph_logfc=",lfc,"qvalue=",pval,"html",sep=".")
  
  #rmarkdown::render('/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/res_04_02_2021/template_OS_Adi.Rmd', output_file = "report_survival_analysis_coxph_02_24_2021.html")
  
  #rmarkdown::render('/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/res_04_02_2021/template_OS_Adi.Rmd', output_file = output_file)
  
  setwd("/mnt/data/lab/gmt_data/data_brain/adi_february_2021")
  
}	
