
heatmapMouse<-function(matrix_to_plot_ms,output_mm,ann_mm_clean){
  

              cal_z_score <- function(x){
                (x - mean(x)) / sd(x)
              }

	      colours<-list("Condition"=c("C1"="red","C2"="green","C1A"="orange","C1B"="blue","C2"="yellow","C3"="brown","C8"="violet"),			    "Condition2"=c("C1"="red","C2"="green","C3"="grey","C8"="violet"))

              colAnn <- HeatmapAnnotation(df=ann_mm_clean, which="col",col=colours,                                         
					  annotation_width=unit(c(1, 4), "cm"), 					   
					  gap=unit(1, "mm"),								
					  annotation_name_side = "left",																	  border=TRUE)

              mat_mm<-log(matrix_to_plot_ms+1,2)
             
	      #cor_matrix<-cor(mat_mm)

	      #png(paste("COR",output_mm,sep="."), height=8,width=10,units='in',res=300)

              #corh<-Heatmap(cor_matrix,	
	      #	      top_annotation=colAnn,	         
	      #	      clustering_distance_rows = "pearson",		            
	      #	      clustering_method_rows = "average",			              
	      #	      clustering_distance_columns = "pearson",				         
	      #	      clustering_method_columns = "average")
	      #draw(corh)
	      #dev.off()

	      
	      #png(paste("One_minusCOR",output_mm,sep="."), height=8,width=10,units='in',res=300)
	      #corh<-Heatmap(1-cor_matrix,       
	      #	      top_annotation=colAnn,             
	      #	      clustering_distance_rows = "pearson",                         
	      #	      clustering_method_rows = "average",                                     
	      #	      clustering_distance_columns = "pearson",                                   
	      #	      clustering_method_columns = "average")
	      #draw(corh)
	      #dev.off()

	      
	      
	      mat_mm<-t(apply(mat_mm,1,cal_z_score))
	      colnames(mat_mm)<-colnames(matrix_to_plot_ms)
	     	 
              library(circlize)	      
      	      col_fun = colorRamp2(c(min(ann_mm_clean$scores_global), 0, max(ann_mm_clean$scores_global)), c("deepskyblue3", "gainsboro", "darkorange2"))			    
	      col_fun2 = colorRamp2(c(min(ann_mm_clean$gsva_scores), 0, max(ann_mm_clean$gsva_scores)), c("deepskyblue3", "gainsboro", "darkorange2"))    

	      colours<-list("Condition"=c("C1"="red","C2"="green","C1A"="orange","C1B"="blue","C2"="yellow","C3"="brown","C8"="violet"),
                            "Condition2"=c("C1"="red","C2"="green","C3"="grey","C8"="violet"),
	                    "scores_global"=col_fun,
			    "gsva_scores"=col_fun2)

              colAnn <- HeatmapAnnotation(df=ann_mm_clean, which="col",col=colours,
                                          annotation_width=unit(c(1, 4), "cm"),
                                          gap=unit(1, "mm"),
                                          annotation_name_side = "left",
                                          border=TRUE)


              png(output_mm, height=12,width=8,units='in',res=300)
              
              mmmap<-Heatmap(mat_mm,
                             show_column_names=T,
                             show_row_names=F,
                             cluster_columns = T,
                             show_column_dend = TRUE,
                             row_dend_reorder = F,
                             column_dend_reorder = F,
                             show_row_dend = TRUE,
                             clustering_distance_rows = "spearman",
                             clustering_method_rows = "complete",
                             clustering_distance_columns = "spearman",
                             clustering_method_columns = "complete",
                             top_annotation=colAnn)
              
              draw(mmmap, heatmap_legend_side="bottom", annotation_legend_side="right")
              
              dev.off()

	                    
	      #cor_matrix<-cor(mat_mm)

	      
	      #png(paste("COR_ZSCORE",output_mm,sep="."), height=8,width=10,units='in',res=300)
	      #corh<-Heatmap(cor_matrix,		           
	      #		    top_annotation=colAnn,				
	      #	  	    clustering_distance_rows = "pearson",			
	      #	    	    clustering_method_rows = "average",				
	      #		    clustering_distance_columns = "pearson",				
	      #		    clustering_method_columns = "average")
	      #draw(corh)
	      #dev.off()
}
