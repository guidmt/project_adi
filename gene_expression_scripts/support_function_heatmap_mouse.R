heatmapMouseFromHSclusters<-function(genes_in_clusters_df_var,output){

	for(clt in unique(genes_in_clusters_df_var[,2])){

	output_mm<-paste('Mus_musculus',output,clt,".png",sep='')
	
	genes_in_clusters_hs<-genes_in_clusters_df_var[which(genes_in_clusters_df_var[,2]%in%clt),1]

	genes_in_clusters_mm<-ortholog_table_res[ortholog_table_res[,1]%in%genes_in_clusters_hs,2]

	matrix_mm_ort_with_gene<-data.frame(genes=rownames(matrix_to_plot_ms_ort),matrix_to_plot_ms_ort)

	heatmapMouse(matrix_mm_ort_with_gene,
		     output_mm=output_mm,
		     subset_genes=TRUE,
		     genes_to_select=genes_in_clusters_mm,
		     row_names=T)
}
}


heatmapMouse<-function(matrix_to_plot_mouse,output_mm,subset_genes=F,genes_to_select,row_names=F){

cal_z_score <- function(x){
	  (x - mean(x)) / sd(x)
}


if(subset_genes==TRUE){

mat_mm<-log(matrix_to_plot_mouse[which(matrix_to_plot_mouse[,1]%in%genes_to_select),-1]+1,2)
mat_mm<-t(apply(mat_mm,1,cal_z_score))

}

if(subset_genes==FALSE){

mat_mm<-log(matrix_to_plot_mouse[,-1]+1,2)
mat_mm<-t(apply(mat_mm,1,cal_z_score))

}

mat_mm[is.na(mat_mm)]<-0

col_fun = colorRamp2(c(min(ann_mm_clean$scores_in_mouse), 0, max(ann_mm_clean$scores_in_mouse)), c("deepskyblue3", "gainsboro", "darkorange2"))

colours<-list("Condition"=c("C1"="red","C2"="green","C1A"="orange","C1B"="blue","C2"="yellow","C3"="brown","C8"="violet"),
	                    "Condition2"=c("C1"="red","C2"="green","C3"="grey","C8"="violet"),
			    col=list(scores_in_mouse=col_fun))

colAnn <- HeatmapAnnotation(df=ann_mm_clean, which="col",col=colours,
                            annotation_width=unit(c(1, 4), "cm"),
                            gap=unit(1, "mm"),
                            annotation_name_side = "left",
                            border=TRUE)

png(output_mm, height=8,width=8,units='in',res=300)

mmmap<-Heatmap(mat_mm,
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

draw(mmmap, heatmap_legend_side="bottom", annotation_legend_side="right")

dev.off()

}
