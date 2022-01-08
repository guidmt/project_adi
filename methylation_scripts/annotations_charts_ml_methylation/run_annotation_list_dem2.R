library(GenomicRanges)

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis")
source("support_function_heatmap_mouse.ge.R")

setwd("/mnt/data/lab/gmt_data/data_brain/adi/BriefSummaryReports")

input_met_mm_rid<-read.delim(file="/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/regions_dmrs_adi_pval0.05_lfc1.txt",header=T)
colnames(input_met_mm_rid)<-c("chr","Cpg_start","Cpg_end","gene_symbol")
input_met_mm_rid[,1]<-gsub(input_met_mm_rid[,1],pattern="chr",replacement="")

df_ann<-unique(input_met_mm_rid)

#
# which are the methylation levels in mouse?
#

library(data.table)

setwd("/mnt/data/lab/gmt_data/data_brain/adi/CytosineReports")

list_meth_file<-dir()

list_samples<-vector(mode='list',length(list_meth_file))

for(i in 1:length(list_meth_file)){

	print("reading data")
	meth_table_current_mouse<-data.frame(fread(file=list_meth_file[i]))
	colnames(meth_table_current_mouse)<-c("chromosome","position","strand","count_methylated","count_unmethylated","c_context","trinucleotide_context")

	beta_current_values<-data.frame()

	#select only the sites in which the gens were identified DEGs
	for(row_ann in 1:nrow(df_ann)){
	
		chrom<-paste("chr",df_ann[row_ann,1],sep="")
		start<-df_ann[row_ann,2]
		end<-df_ann[row_ann,3]
		genes<-df_ann[row_ann,ncol(df_ann)]
	

		test<-meth_table_current_mouse[which(meth_table_current_mouse[,1]==chrom &
    						       meth_table_current_mouse[,2]>=start & meth_table_current_mouse[,2]<=end),]
		
		if(nrow(test)>=1){

		print("true")

		coordinates<-paste(range(test[,1]),collapse="-")

		#https://www.biostars.org/p/445286/
		#https://github.com/FelixKrueger/Bismark/issues/354
		#https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-10-r83

		beta_values<-sum(test$count_methylated)/(sum(test$count_methylated)+sum(test$count_unmethylated))

		df_beta<-data.frame(gene_symbol=genes,beta_values)
		
		beta_current_values<-rbind(beta_current_values,df_beta)

		}
	
	}


list_samples[[i]]<-beta_current_values

}


res_meth<-do.call(cbind,list_samples)

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

library("ComplexHeatmap")

png("mm_heatmap_regions_dmrs_adi_pval0.05_lfc1.newestimation.png", height=8,width=8,units='in',res=300)

hsmatrix<-as.matrix(res_meth3[,-1])
fx<-function(x){(x-mean(x))/(max(x)-min(x))}
hsmatrix<-t(apply(hsmatrix,1,fx))

hsmap<-Heatmap(hsmatrix,
               show_column_names=TRUE,
               show_row_names=FALSE,
               cluster_columns = T,
               show_column_dend = TRUE,
               row_dend_reorder = FALSE,
               column_dend_reorder = FALSE, 
	       show_row_dend = TRUE,		
	       clustering_distance_columns = "euclidean",
	       clustering_method_columns = "complete",
	       col=circlize::colorRamp2(c(0,0.5,1), c("Darkblue", "orange", "red")))
draw(hsmap)

dev.off()

#
# Plot methylation without a sample that is repeat and no control
#

png("mm_heatmap_regions_dmrs_adi_pval0.05_lfc1.newestimation_no_redo_no_control.png", height=8,width=10,units='in',res=300)

columns_to_save<-grep(grep(grep(colnames(res_meth3),pattern="redo",invert=T,value=T),pattern="Control",invert=T,value=T),pattern="gene_symbol",invert=T,value=T)

ha = HeatmapAnnotation(
		       status = c("MUT","MUT","MUT","IWT","IWT","IWT"),
	               col = list(
				  status = c("MUT" = "red", 
					     "MUT" = "red", 
					     "MUT" = "red",
					     "IWT"  = "green2",
					     "IWT"  = "green2",
					     "IWT"  = "green2")
				  )			   
		       )

hsmatrix<-as.matrix(res_meth3[,columns_to_save])
fx<-function(x){(x-mean(x))/(max(x)-min(x))}
hsmatrix<-t(apply(hsmatrix,1,fx))
hsmatrix[is.na(hsmatrix)]<-0 #there are few values with NA

hsmap<-Heatmap(hsmatrix,
               show_column_names=TRUE,
               show_row_names=FALSE,
               cluster_columns = T,
               show_column_dend = TRUE,
               row_dend_reorder = FALSE,
               column_dend_reorder = FALSE,
               show_row_dend = TRUE,
               clustering_distance_columns = "euclidean",
	       clustering_method_columns = "complete",
	       top_annotation=ha,
	       column_split=factor(c("MUT","MUT","MUT","IWT","IWT","IWT"),levels=c("MUT","IWT")),	         
     	       col=circlize::colorRamp2(c(0,0.5,1), c("Darkblue", "orange", "red"))
	       )

draw(hsmap)

dev.off()


#
# Plot methylation without a sample that is repeat and control, split by groups
#

png("mm_heatmap_regions_dmrs_adi_pval0.05_lfc1.newestimation_no_redo_yes_control.png", height=8,width=10,units='in',res=300)

columns_to_save<-grep(grep(colnames(res_meth3),pattern="redo",invert=T,value=T),pattern="gene_symbol",invert=T,value=T)

status<-factor(c("CTRL","MUT","MUT","MUT","IWT","IWT","IWT"),levels=c("IWT","CTRL","MUT"))

ha = HeatmapAnnotation(
                       status = status,
                       col = list(
                       status = c("CTRL"="blue",
				  "MUT" = "red", 
                                  "MUT" = "red", 
                                  "MUT" = "red",
                                  "IWT"  = "green2",
				  "IWT"  = "green2",
				  "IWT"  = "green2")
		       )                        
		       )

hsmatrix<-as.matrix(res_meth3[,columns_to_save])

fx<-function(x){(x-mean(x))/(max(x)-min(x))}
hsmatrix<-t(apply(hsmatrix,1,fx))

hsmap<-Heatmap(hsmatrix,
                         show_column_names=TRUE,
			 show_row_names=FALSE,
			 show_column_dend = TRUE,
			 row_dend_reorder = FALSE,
			 show_row_dend = TRUE,
			 clustering_distance_columns = "euclidean",
			 clustering_distance_rows = "euclidean",			 
			 clustering_method_rows = "average",
			 clustering_method_columns = "average",
			 top_annotation=ha,
			 column_split=status,               
			 col=circlize::colorRamp2(c(0,0.5,1), c("Darkblue", "orange", "red"))
			 )
hsmap

dev.off()

save(res_meth3,file="methylation_data_DMET_CpG_05_11_2020_pval0.05_lfc1_newestimation.RData")

run_after=F

if(run_after==TRUE){


#
#Are these genes hyper o hypo methylated?
#

check_iper_ipo<-res_meth3[,-1]

#
# IAP89, IAP89B, IAPA and IwtAP2R4a
#
ratio_beta<-apply(check_iper_ipo[,c(2:4)],1,mean)-apply(check_iper_ipo[,c(5,7,8)],1,mean)
ratios<-ifelse(ratio_beta>0.2,"iper","ipo")

dfres<-data.frame(genes=res_meth3[,1],
                  ratio_beta=ratio_beta,
                  status=ratios,
                  mu_mut=apply(check_iper_ipo[,c(2:4)],1,mean),
                  mu_wt=apply(check_iper_ipo[,c(5,7,8)],1,mean),
                  ratio_m=ratio_beta)

write.table(dfres,"methylation_data_DMET_CpG.STATUS.2000_10_07_2020.noredo.nocontrol.newestimation.txt",sep="\t",row.names=F)

iper_genes<-dfres[dfres$status=="iper","genes"]

png("mm_meth_degs.1000.newestimation_no_redo_no_control.only_iper_genes.png", height=8,width=10,units='in',res=300)

columns_to_save<-grep(grep(grep(colnames(res_meth3),pattern="redo",invert=T,value=T),pattern="Control",invert=T,value=T),pattern="gene_symbol",invert=T,value=T)

ha = HeatmapAnnotation(
  status = c("MUT","MUT","MUT","IWT","IWT","IWT"),
  col = list(
    status = c("MUT" = "red",
               "MUT" = "red",
               "MUT" = "red",
               "IWT"  = "green2",
               "IWT"  = "green2",
               "IWT"  = "green2")
  )
)

hsmap<-Heatmap(as.matrix(res_meth3[iper_genes,columns_to_save]),
               show_column_names=TRUE,
               show_row_names=FALSE,
               cluster_columns = T,
               show_column_dend = TRUE,
               row_dend_reorder = FALSE,
               column_dend_reorder = FALSE,
               show_row_dend = TRUE,
               clustering_distance_columns = "euclidean",
               clustering_method_columns = "complete",
               top_annotation=ha,
               column_split=factor(c("MUT","MUT","MUT","IWT","IWT","IWT"),levels=c("MUT","IWT"))
               
)

draw(hsmap)

dev.off()


#
# what about the impact of the methylation levels on the the gene-expression levels of the genes?
#

# read sampleinfo file
setwd("/mnt/data/lab/gmt_data/data_brain/")
mm_sampleinfo<-read.delim(file='sampleinfo2.txt')

mm_sampleinfo2<-mm_sampleinfo[which(mm_sampleinfo$Condition2%in%c("C1","C2")),]
mm_sampleinfo2<-mm_sampleinfo2[grep(mm_sampleinfo2$Run,pattern="CSHL",invert=T),]

# read the expression data 
setwd("/mnt/data/lab/gmt_data/data_brain/adi_26_03_2020/rawdata")
tab_expression_mouse<-read.csv(file="Adi_IDH.tpm.csv")
tab_expression_mouse<-tab_expression_mouse[which(tab_expression_mouse$Gene_name%in%res_meth3[,1]),]

tab_exp_mm<-data.frame(gs=tab_expression_mouse$Gene_name,
		       tab_expression_mouse[,which(colnames(tab_expression_mouse)%in%mm_sampleinfo2$File_name)])

tab_exp_mm2<-aggregate(tab_exp_mm,by=list(tab_exp_mm$gs),FUN=mean)[,-2]
rownames(tab_exp_mm2)<-tab_exp_mm2[,1]

setwd("/mnt/data/lab/gmt_data/data_brain/")
mm_sampleinfo<-read.delim(file='sampleinfo2.txt')
mm_sampleinfo_rid<-mm_sampleinfo[mm_sampleinfo$Condition %in% c("C1A","C1B","C2"),]
mm_sampleinfo2<-mm_sampleinfo_rid[,c('File_name','Condition','Condition2')]
rownames(mm_sampleinfo2)<-mm_sampleinfo2[,1]

tab_exp_mm2<-tab_exp_mm2[,-1]
tab_exp_mm3<-tab_exp_mm2[,colnames(tab_exp_mm2)%in% mm_sampleinfo_rid$File_name]

mm_sampleinfo2<-mm_sampleinfo2[mm_sampleinfo2$File_name%in%colnames(tab_exp_mm2),]
mm_sampleinfo2<-mm_sampleinfo2[match(colnames(tab_exp_mm2),rownames(mm_sampleinfo2)),]

output_mm<-"GE_CpG_diff_mut_vs_wt_in_expressionC1_C2.noredo.nocontrol.newestimation.1000_10_07_2020.png"

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis")

png(output_mm, height=8,width=10,units='in',res=300)

cal_z_score <- function(x){
          (x - mean(x)) / sd(x)
}

mat_mm<-log(tab_exp_mm3+1,2)
mat_mm<-t(apply(mat_mm,1,cal_z_score))

mat_mm[is.na(mat_mm)]<-0

colours<-list("Condition"=c("C1A"="orange","C1B"="blue","C2"="yellow"),
	                                  "Condition2"=c("C1"="red","C2"="green"))

colAnn <- HeatmapAnnotation(df=mm_sampleinfo2[,-1], which="col",col=colours,
                            annotation_width=unit(c(1, 4), "cm"),
                            gap=unit(1, "mm"),
                            annotation_name_side = "left",
		            border=TRUE)


mmmap<-Heatmap(mat_mm,
	       show_row_names=F,
	       show_column_names=FALSE,
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

}
