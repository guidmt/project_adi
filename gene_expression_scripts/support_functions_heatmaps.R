plot_heatmap_hs<-function(matrix_to_plot_hs,annotation,output_hs,output_hs2,column_title = "Title",row_names=F,split=F,how_split=NULL,width=8,height=7,pos_legend_heatmap="bottom",pos_legend_ann="right",km=FALSE,number_clusters=3,number_clusters_row=3,merged_with_array=FALSE,nopval=FALSE,strong_cluster=FALSE,png_output=T){

  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
 

library(circlize)
col_fun = colorRamp2(c(min(annotation$scores_global), 0, max(annotation$scores_global)), c("deepskyblue3", "gainsboro", "darkorange2"))
col_fun2 = colorRamp2(c(min(annotation$gsva_scores), 0, max(annotation$gsva_scores)), c("deepskyblue3", "gainsboro", "darkorange2"))

colours <- list("type"=c("GBM"="red2",
	 "LGG"="skyblue1",
	 "no_ann"="gainsboro",
	 "Brain - Cortex"="violet",
	 "Brain - Cerebellum"="orange",
	 "Brain - Frontal Cortex (BA9)"="green",
	 "Brain - Caudate (basal ganglia)"="cyan",
	 "Brain - Nucleus accumbens (basal ganglia)"="yellow2",
	 "Brain - Putamen (basal ganglia)"="blue",
	 "Brain - Hypothalamus"="gold",
	 "Brain - Hypothalamus"="cornflowerblue",
	 "Brain - Spinal cord (cervical c-1)"="tan",
	 "Brain - Hippocampus"="snow4",
	 "Brain - Anterior cingulate cortex (BA24)"="sienna1",
	 "Brain - Cerebellar Hemisphere"="plum1",
	 "Brain - Substantia nigra"="dimgrey",
	 "Brain - Amygdala"="violetred3"
	 ),
	  "clusters"=c("LGr1"="brown2","LGr2"="dodgerblue4","LGr3"="darkorchid3","LGr4"="olivedrab3","no_class"="grey","no_ann"="white"),
	  "ATRX"=c("1"="royalblue1","0"="gainsboro","no_ann"="white"),
	  "IDH1"=c("1"="red","0"="gainsboro","no_ann"="white"),
	  "IDH1_R132H"=c("1"="red","0"="gainsboro","no_ann"="white"),
	  "IDH1_other"=c("1"="yellow","0"="gainsboro","no_ann"="white"),
	  "PDGFRA"=c("1"="green","0"="gainsboro","no_ann"="white"),
	  "TP53"=c("1"="violet","0"="gainsboro","no_ann"="white"),
	  "histology"=c("astrocytoma"="slateblue4",
		  "glioblastoma"="violetred2",
		  "no_histology"="khaki3",
		  "oligoastrocytoma"="limegreen",
		  "oligodendroglioma"="goldenrod2",
		  "no_ann"="gainsboro",
		  "no_histology"="gainsboro"),
		"pdgfra_exp"=c("up"="firebrick2","down"="dodgerblue1","no_ann"="gainsboro"),
		"atrx_exp"=c("up"="gold","down"="darkslategray4","no_ann"="gainsboro"),
		"1p_19q_deletion"=c("codel"="cadetblue4","non-codel"="gainsboro","no_ann"="white"),
		"IDH_codel_subtype"=c("IDHmut-codel"="darkorange3","IDHmut-non-codel"="gold2","IDHwt"="chartreuse4","no_ann"="white"),
		"chr7_gain_chr_10loss"=c("Gain_chr_7_&_loss_chr_10"="darkslateblue","No_combined_CNA"="gainsboro","no_ann"="white"),
		"chr19_20_cogain"=c("Gain_chr_19/20"="deeppink2","No_chr_19/20_gain"="gainsboro","no_ann"="white"),
		"TERT_promoter_status"=c("Mutant"="turquoise1","WT"="gainsboro","no_ann"="white"),
	"experimental_model"=c("C1"="firebrick1","C2"="dodgerblue1","C3"="yellow3","C8"="darkorange1","C4"="forestgreen","no_ann"="gainsboro","other"="gainsboro"),
	"IDH1_status_microarray"=c("1"="red","0"="yellow","no_ann"="gainsboro"),
	"Grade"=c("WHO IV"="red2","WHO III"="orange2","WHO II"="green2","no_ann"="gainsboro","WHO I"="blue"),
	"1p_19q_codeletion"=c("Codel"="violet","Non-codel"="blue2","no_ann"="gainsboro"),
	"study"=c("TCGA"="blue2","CGGA"="red2","ann_GSE108474"="pink","ann_GSE16011"="green2","ann_GSE70231"="orange2"),
	"TME"=c("Cytotoxic/Exhausted"="red2","MixedImmunity"="deepskyblue2","ImmuneDesert"="burlywood3"),
	"status_C1AC1B"=c("high"="orange2","low"="blue2"),
	"mod2_classification"=c("cl5"="red2","cl4"="orange2","cl3"="green2","cl2"="violetred2","cl1"="slateblue4"),
	"scores_global"=col_fun,
	"gsva_scores"=col_fun2)

#sample_for_analysis<-annotation[,1]

#select only the samples for the current analysis
#matrix_to_plot_hs<-matrix_to_plot_hs[,which(colnames(matrix_to_plot_hs)%in%sample_for_analysis)]

#matrix_to_plot_hs<-matrix_to_plot_hs[,match(colnames(matrix_to_plot_hs),annotation[,1])]
#match the samples between colnames of the matrix and annotation file

mat_temp_trs<-t(matrix_to_plot_hs)
mat_temp_trs_with_id<-data.frame(ID=rownames(mat_temp_trs),mat_temp_trs)

size_ann<-ncol(annotation)
  
merged_exp_ann<-merge(annotation,mat_temp_trs_with_id,by.x="Samples",by.y="ID")

annotation<-merged_exp_ann[,c(1:size_ann)]
matrix_to_plot_hs<-t(merged_exp_ann[,-c(1:size_ann)])
colnames(matrix_to_plot_hs)<-merged_exp_ann$Samples


print(dim(annotation))
print(dim(matrix_to_plot_hs))

colAnn <- HeatmapAnnotation(df=annotation[,-1], which="col",col=colours,
			    annotation_width=unit(c(1, 4), "cm"),
			    gap=unit(1, "mm"),
			    annotation_name_side = "left",
			    border=TRUE)

#
# HS: plot the heatmap of the genes found in mm in hs 
#

if(png_output==TRUE){
  
png(output_hs, height=8,width=width,units='in',res=300)

}else{

  pdf(output_hs,width=width,height=height)
}

if(merged_with_array==FALSE){

mat_hs<-log(matrix_to_plot_hs+1,2)
} else{

mat_hs<-matrix_to_plot_hs

}
output_tummap<-paste(sapply(strsplit(output_hs,split=".png"),"[[",1),"tumap.txt",sep=".")
write.table(mat_hs,file=output_tummap,sep="\t",row.names=F,quote=F)

mat_hs<-t(apply(mat_hs,1,cal_z_score))
colnames(mat_hs)<-colnames(matrix_to_plot_hs)

#
# If the expression of a gene is 0 on TCGA can return NA after z-score
#
mat_hs<-na.omit(mat_hs)
print("Size unfiltered heatmap")
print(dim(mat_hs))

if(split==FALSE){
  
hsmap<-Heatmap(mat_hs,
            	 show_column_names=FALSE,
            	 show_row_names=row_names,
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

draw(hsmap, heatmap_legend_side=pos_legend_heatmap, annotation_legend_side=pos_legend_ann,column_title = column_title)

}


if(split==TRUE){

	if(km!=TRUE){

	hsmap<-Heatmap(mat_hs,
      		       show_column_names=FALSE,
      	               show_row_names=row_names,
      		       cluster_columns = T,
      		       show_column_dend = TRUE,
      		       row_dend_reorder = F,
      		       column_dend_reorder = F,
      		       show_row_dend = TRUE,
      		       clustering_distance_rows = "euclidean",
      		       clustering_method_rows = "complete",
      		       clustering_distance_columns = "euclidean",
      		       clustering_method_columns = "complete",
      		       top_annotation=colAnn,
		       column_split=annotation[,how_split],
		       cluster_column_slices = FALSE,
		       row_km=number_clusters_row)	
	
	} else {
	
	hsmap<-Heatmap(mat_hs,
      	         show_column_names=FALSE,
      	         show_row_names=row_names,
      	         cluster_columns = T,
		 show_column_dend = TRUE,
      		 row_dend_reorder = F,
      		 column_dend_reorder = F,
      		 show_row_dend = TRUE,
      		 clustering_distance_rows = "euclidean",
      		 clustering_method_rows = "complete",
      	         clustering_distance_columns = "euclidean",
      	         clustering_method_columns = "complete",
      	         top_annotation=colAnn,
      	         column_km=number_clusters,
      		 row_km=number_clusters_row)
	}
  
	hm<-draw(hsmap, heatmap_legend_side=pos_legend_heatmap, annotation_legend_side=pos_legend_ann,column_title = column_title)
  	print(hm)
	
       #get the orders of the genes

    if(km==TRUE | km==FALSE & strong_cluster==TRUE){
      
       for (i in 1:length(row_order(hm))){

       if (i == 1) {
             clu <- t(t(row.names(mat_hs[row_order(hm)[[i]],])))
             out <- cbind(clu, paste("cluster",  names(row_order(hm)[i]), sep=""))
             colnames(out) <- c("GeneID", "Cluster")
      
       } else {

             clu <- t(t(row.names(mat_hs[row_order(hm)[[i]],])))
      
             clu <- cbind(clu, paste("cluster", names(row_order(hm)[i]), sep=""))
             out <- rbind(out, clu)
       }


       }

       print("check clustering  no filter")
       print("second check")
       print(length(unlist(row_order(hm))))

       print(dim(out))
       print(length(unique(out[,1])))

       write.table(file=paste(output_hs,".clusteringRows.txt",sep=""),out,sep="\t")
  
}

dev.off()

}

#
# HS: plot the heatmap of the genes found in mm in hs, filtering the low variances genes 
#

#variance_genes<-rowVars(matrix_to_plot_hs)

#variance_genes_sorted<-variance_genes[order(variance_genes)]

#cut_off<-quantile(variance_genes_sorted)[2]

#genes_after_variance_filter<-names(variance_genes_sorted[variance_genes_sorted>=cut_off])

#matrix_to_plot_hs_var<-matrix_to_plot_hs[which(rownames(matrix_to_plot_hs)%in%genes_after_variance_filter),]

#if(merged_with_array==FALSE){

#mat_hs_var<-log(matrix_to_plot_hs_var+1,2)

#} else {

#mat_hs_var<-matrix_to_plot_hs_var

#}	

#output_tummap<-paste(sapply(strsplit(output_hs2,split=".png"),"[[",1),"tumap.txt",sep=".")
#write.table(mat_hs_var,file=output_tummap,sep="\t",row.names=F,quote=F)

#mat_hs_var<-t(apply(mat_hs_var,1,cal_z_score))
#colnames(mat_hs_var)<-colnames(matrix_to_plot_hs)

#mat_hs_var<-na.omit(mat_hs_var)

#print("Size filtered heatmap")
#print(dim(mat_hs_var))

#png(output_hs2, height=8,width=width,units='in',res=300)

#if(split==FALSE){

#print(split)
  
#hsmap<-Heatmap(mat_hs_var,
#            	 show_column_names=FALSE,
#          	   show_row_names=row_names,
#            	 cluster_columns = T,
#          	   show_column_dend = TRUE,
#          	   row_dend_reorder = F,
#          	   column_dend_reorder = F,									       
#            	 show_row_dend = TRUE,
#          	   clustering_distance_rows = "euclidean",
#          	   clustering_distance_columns = "euclidean",
#          	   clustering_method_columns = "complete",
#          	   top_annotation=colAnn)
#
#draw(hsmap, heatmap_legend_side=pos_legend_heatmap, annotation_legend_side=pos_legend_ann,column_title = column_title)

#}


#if(split==TRUE){

#if(km!=TRUE){
	
#hsmap<-Heatmap(mat_hs_var,
#      	       show_column_names=FALSE,
#      	       show_row_names=row_names,
#      	       cluster_columns = T,
#      	       show_column_dend = TRUE,
#      	       row_dend_reorder = F,
#      	       column_dend_reorder = F,
#      	       show_row_dend = TRUE,
#      	       clustering_distance_rows = "euclidean",
#      	       clustering_method_rows = "complete",
#      	       clustering_distance_columns = "euclidean",
#      	       clustering_method_columns = "complete",
#      	       top_annotation=colAnn,
#      	       column_split=how_split)
      
#}else{
	 
 	
#hsmap<-Heatmap(mat_hs_var,
#      	       show_column_names=FALSE,
#      	       show_row_names=row_names,
#      	       cluster_columns = T,
#      	       show_column_dend = TRUE,
#      	       row_dend_reorder = F,
#      	       column_dend_reorder = F,
#      	       show_row_dend = TRUE,
#      	       clustering_distance_rows = "euclidean",
#      	       clustering_method_rows = "complete",
#      	       clustering_distance_columns = "euclidean",
#      	       clustering_method_columns = "complete",
#      	       top_annotation=colAnn,
#      	       column_km=number_clusters,
#      	       row_km=number_clusters_row)
#}
  
#  hmvar<-draw(hsmap, heatmap_legend_side=pos_legend_heatmap, annotation_legend_side=pos_legend_ann,column_title = column_title)
  
#  print(hmvar)
  
  #
  # Get the order of the genes
  #
#  if(km==TRUE){
    
#  out<-data.frame()
       
#  for (i in 1:length(row_order(hmvar))){

#  if (i == 1) {
		        
#  clu <- t(t(row.names(mat_hs[row_order(hmvar)[[i]],])))
#  out <- cbind(clu, paste("cluster", names(row_order(hmvar)[i]), sep=""))
#  colnames(out) <- c("GeneID", "Cluster")

#  } else {

#  clu <- t(t(row.names(mat_hs[row_order(hmvar)[[i]],])))

#  clu <- cbind(clu, paste("cluster", names(row_order(hmvar)[i]), sep=""))
#  out <- rbind(out, clu)
#       	}
	 
#  }
       
#       write.table(file=paste(output_hs2,".clusteringRows.txt",sep=""),out,sep="\t")
       
       
#  }
  
#}

#dev.off()

}

