plot_heatmap_hs<-function(matrix_to_plot_hs,annotation,output_hs,output_hs2,column_title = "Title"){

  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  
colours <- list("type"=c("GBM"="red2","LGG"="green2"),
                "clusters"=c("LGr1"="brown2","LGr2"="dodgerblue4","LGr3"="darkorchid3","LGr4"="olivedrab3","no_class"="grey"),
                "ATRX"=c("1"="royalblue1","0"="gainsboro"),
                "IDH1_R132H"=c("1"="red","0"="gainsboro"),
                "IDH1_other"=c("1"="yellow","0"="gainsboro"),
                "PDGFRA"=c("1"="green","0"="gainsboro"),
                "TP53"=c("1"="violet","0"="gainsboro"),
                "histology"=c("astrocytoma"="slateblue4",
                              "glioblastoma"="violetred2",
                              "no_histology"="khaki3",
                              "oligoastrocytoma"="limegreen",
                              "oligodendroglioma"="goldenrod2")
)

colAnn <- HeatmapAnnotation(df=annotation[,-1], which="col",col=colours,
                            annotation_width=unit(c(1, 4), "cm"),
                            gap=unit(1, "mm"),
                            annotation_name_side = "left",
                            border=TRUE)

sample_for_analysis<-annotation[,1]

#select only the samples for the current analysis
matrix_to_plot_hs<-matrix_to_plot_hs[,which(colnames(matrix_to_plot_hs)%in%sample_for_analysis)]

matrix_to_plot_hs<-matrix_to_plot_hs[,match(colnames(matrix_to_plot_hs),annotation[,1])]
#match the samples between colnames of the matrix and annotation file


#
# HS: plot the heatmap of the genes found in mm in hs 
#

png(output_hs, height=8,width=8,units='in',res=300)

mat_hs<-log(matrix_to_plot_hs+1,2)

mat_hs<-t(apply(mat_hs,1,cal_z_score))
colnames(mat_hs)<-colnames(matrix_to_plot_hs)

hsmap<-Heatmap(mat_hs,
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

draw(hsmap, heatmap_legend_side="bottom", annotation_legend_side="right",column_title = column_title)

dev.off()


#
# HS: plot the heatmap of the genes found in mm in hs, filtering the low variances genes 
#

variance_genes<-rowVars(matrix_to_plot_hs)

variance_genes_sorted<-variance_genes[order(variance_genes)]

cut_off<-quantile(variance_genes_sorted)[2]

genes_after_variance_filter<-names(variance_genes_sorted[variance_genes_sorted>=cut_off])

matrix_to_plot_hs_var<-matrix_to_plot_hs[which(rownames(matrix_to_plot_hs)%in%genes_after_variance_filter),]

mat_hs_var<-log(matrix_to_plot_hs_var+1,2)

mat_hs_var<-t(apply(mat_hs_var,1,cal_z_score))
colnames(mat_hs_var)<-colnames(matrix_to_plot_hs)

png(output_hs2, height=8,width=8,units='in',res=300)

hsmap<-Heatmap(mat_hs_var,
               show_column_names=FALSE,
               show_row_names=F,
               cluster_columns = T,
               show_column_dend = TRUE,
               row_dend_reorder = F,
               column_dend_reorder = F,                                                                                                                                
               show_row_dend = TRUE,
               clustering_distance_rows = "euclidean",
               clustering_distance_columns = "euclidean",
               clustering_method_columns = "complete",
               top_annotation=colAnn)
draw(hsmap, heatmap_legend_side="bottom", annotation_legend_side="right",column_title= paste(column_title,"variance filter"))

dev.off()

}
