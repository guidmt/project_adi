plot_heatmap_hs_ge_meth<-function(matrix_to_plot_hs,mat_meth,annotation,output_hs,output_hs2,height=8,output_hs3,qthr=10,column_title = "Title",row_names=F,split=T,how_split=NULL,width=8,pos_legend_heatmap="bottom",pos_legend_ann="right",km=FALSE,number_clusters=3,number_clusters_row=3,mut_trim=200,variables_to_use_tsne=c("type","IDH1_R132H","IDH_methylation_clusters","histology"),mut_burden=F,perplex=15,tsne=F,number_clusters_row_meth=3){

library(circlize)
library(stringr)
library(stringi)

current_dir<-getwd()

  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  
# Upload the colors
source("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/colour_annotation.R")
sample_for_analysis<-annotation[,1]

#select only the samples for the current analysis
matrix_to_plot_hs<-matrix_to_plot_hs[,which(colnames(matrix_to_plot_hs)%in%sample_for_analysis)]

#matrix_to_plot_hs<-matrix_to_plot_hs[,match(colnames(matrix_to_plot_hs),annotation[,1])]

#
# HS: plot the heatmap of the genes found in mm in hs 
#

mat_hs<-log(matrix_to_plot_hs+1,2)

# mat_hs_for_tsne<-log(matrix_to_plot_hs+1,2)

# output_tummap<-paste(sapply(strsplit(output_hs,split=".png"),"[[",1),"tumap.txt",sep=".")
# write.table(mat_hs,file=output_tummap,sep="\t",row.names=F,quote=F)

mat_hs<-t(apply(mat_hs,1,cal_z_score))
colnames(mat_hs)<-colnames(matrix_to_plot_hs)

#
# get the common samples between gene-expression adn methylation data.
#

samples_ge<-colnames(mat_hs)
samples_meth<-colnames(mat_meth)[-c(1)]

common_samples<-intersect(intersect(samples_ge,samples_meth),annotation[,1])

print(length(common_samples))

if(split==TRUE & km == FALSE){
print("get samples for split")

how_split<-how_split[names(how_split)%in%common_samples]
print(length(how_split))

}

#
# Define the gene-expression data 
#
mat_hs<-na.omit(mat_hs[,which(colnames(mat_hs)%in%common_samples)])

# colnames(mat_hs)<-gsub(make.names(colnames(mat_hs),unique=T),pattern="\\.",replacement="-")
# 
# colnames(mat_meth)[-c(1:2)]<-gsub(make.names(colnames(mat_meth)[-c(1:2)],unique=T),pattern="\\.",replacement="-")

#
# Define the methylation data 
#

expr_col_fun = colorRamp2(c(0,0.5,1), c("blue", "yellow2", "red"))

mat_meth[is.na(mat_meth)]<-0

print(dim(mat_meth))

colnames(mat_meth)[1]<-"genes"

#
# Select the most variable islands
#

varCpGislands<-apply(mat_meth[,-c(1)],1,var)

mat_meth$varCpGislands <- varCpGislands

mat_meth_temp<-mat_meth[order(mat_meth$varCpGislands,decreasing=T),][1:round(((nrow(mat_meth)*qthr)/100)),]

idx_remove<-which(colnames(mat_meth_temp)%in%"varCpGislands")

mat_meth2<-mat_meth_temp[,-idx_remove]

mat_meth3<-mat_meth2[,-1]
mat_meth3<-mat_meth3[,which(colnames(mat_meth3)%in%common_samples)]

mat_meth3<-mat_meth3[,match(colnames(mat_hs),colnames(mat_meth3))]


#
# Now define the top annotation 
#
annotation2<-annotation[annotation[,1]%in%common_samples,]
annotation2[,1]<-gsub(make.names(annotation2[,1],unique=T),pattern="\\.",replacement="-")
colnames(mat_meth3)<-gsub(colnames(mat_meth3),pattern="\\.",replacement="-")
colnames(mat_hs)<-gsub(colnames(mat_hs),pattern="\\.",replacement="-")

mat_meth3<-mat_meth3[,match(colnames(mat_meth3),annotation2[,1])]
mat_hs<-mat_hs[,match(colnames(mat_hs),annotation2[,1])]
  
colAnn <- HeatmapAnnotation(df=annotation2[,-1], which="col",col=colours,
                             annotation_width=unit(c(1, 4), "cm"),
                             gap=unit(1, "mm"),
                             annotation_name_side = "left",
                             border=TRUE)

print("size exp")
print(dim(mat_hs))
print("size meth")
print(dim(mat_meth3))
print("size ann")
print(dim(annotation2))

# ############################
# # tSNE:
# # 1) get the common samples
# # 2) scale the data between 0-1
# ############################
# 
# scale_max<-function(x){x/max(x)}

#mat_hs_for_tsne2<-mat_hs_for_tsne[,which(colnames(mat_hs_for_tsne)%in%common_samples)]
#mat_hs_for_tsne3<-t(apply(mat_hs_for_tsne2,1,scale_max))

setwd(current_dir)

png(output_hs, height=height,width=width,units='in',res=300)

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
      		       column_split=how_split,
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
  
}	

	if(km==TRUE){

	hsmet<-Heatmap(mat_meth3,
               show_column_names=FALSE,
               show_row_names=row_names,
               show_column_dend = TRUE,
               show_row_dend = TRUE,
               clustering_distance_rows = "euclidean",
               clustering_method_rows = "complete",
               clustering_distance_columns = "euclidean",
               clustering_method_columns = "complete",
               col = expr_col_fun,
	       column_km=number_clusters,
	       row_km=number_clusters_row)

	}else{
	
	hsmet<-Heatmap(mat_meth3,
		       show_column_names=FALSE,
		       show_row_names=row_names,
		       column_dend_reorder = F,
		       show_row_dend = TRUE,
		       clustering_distance_rows = "euclidean",
		       clustering_method_rows = "complete",
		       clustering_distance_columns = "euclidean",
		       clustering_method_columns = "complete",
		       col = expr_col_fun,
		       row_km=number_clusters_row_meth)
		       
	}

	h<-draw(hsmap%v%hsmet)
	
	print(h)

	dev.off()


colnames(mat_hs)<-gsub(make.names(colnames(mat_hs),unique=T),pattern="\\.",replacement="-")
colnames(mat_meth3)<-gsub(colnames(mat_meth3),pattern="\\.",replacement="-")

print(colnames(mat_hs)[1:20])
print(colnames(mat_meth3)[1:20])

mat_tummap<-rbind(mat_hs,mat_meth3)

output_tummap2<-paste(sapply(strsplit(output_hs,split=".png"),"[[",1),"tumap.ge.meth.txt",sep=".")
output_ann_tummap2<-paste(sapply(strsplit(output_hs,split=".png"),"[[",1),"tumap.ge.meth.annotation.txt",sep=".")
# header_to_use<-sapply(strsplit(colnames(mat_tummap),split="-"),FUN=function(X){paste(X[1:3],collapse="-")})

write.table(mat_tummap,file=output_tummap2,sep="\t",row.names=F,quote=F)

write.table(annotation2,file=output_ann_tummap2,sep="\t",row.names=F,quote=F,col.names=T)

if(km==TRUE){

      get_index_columns<-unlist(column_order(h))
      clusters_columns<- rep(names(column_order(h)),unlist(lapply(column_order(h),length)))
      df_columns<-data.frame(index_columns=get_index_columns,columns_names=colnames(mat_hs)[get_index_columns],clusters_columns=clusters_columns)


output_columns2<-paste(sapply(strsplit(output_hs,split=".png"),"[[",1),".clusteringcolumns.ge.meth.txt",sep=".")
write.table(df_columns,file=output_columns2,sep="\t",row.names=F,quote=F)

}

        if(km==TRUE){

        hsmet2<-Heatmap(mat_meth3,
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
                       row_km=number_clusters_row,
		       col=expr_col_fun)
	        
		    png(output_hs2, height=8,width=width,units='in',res=300)	
		    hmetspec<-draw(hsmet2)
	            print(hmetspec)
		    dev.off()

        }

#
# heatmap of only the methylation samples that are in the matrix to complete
#
  
  # mat_meth_all<-mat_meth
  # 
  # varCpGislands<-apply(mat_meth_all[,-c(1:2)],1,var)
  # 
  # mat_meth_all$varCpGislands <- varCpGislands
  # 
  # mat_meth_temp<-mat_meth_all[order(mat_meth_all$varCpGislands,decreasing=T),][1:round(((nrow(mat_meth_all)*qthr)/100)),]
  # 
  # idx_remove<-which(colnames(mat_meth_temp)%in%"varCpGislands")
  # 
  # mat_meth2_all<-mat_meth_temp[,-idx_remove]
  #   
  # mat_meth3_all<-mat_meth2_all[,-c(1:2)]
  # 
  # mat_meth3_all<-mat_meth3_all[,which(colnames(mat_meth3_all)%in%annotation[,1])]
  # 
  # if(km==TRUE){
  #   
  #   colAnn2 <- HeatmapAnnotation(df=annotation[,-1], which="col",col=colours,
  #                                annotation_width=unit(c(1, 4), "cm"),
  #                                gap=unit(1, "mm"),
  #                                annotation_name_side = "left",
  #                                border=TRUE)
  #   
  #   hsmet2<-Heatmap(mat_meth3_all,
  #                   show_column_names=FALSE,
  #                   show_row_names=row_names,
  #                   cluster_columns = T,
  #                   show_column_dend = TRUE,
  #                   row_dend_reorder = F,
  #                   column_dend_reorder = F,
  #                   show_row_dend = TRUE,
  #                   clustering_distance_rows = "euclidean",
  #                   clustering_method_rows = "complete",
  #                   clustering_distance_columns = "euclidean",
  #                   clustering_method_columns = "complete",
  #                   top_annotation=colAnn2,
  #                   column_km=number_clusters,
  #                   row_km=number_clusters_row,
  #                   col=expr_col_fun)
  #   
  #   png(output_hs3, height=8,width=width,units='in',res=300)	
  #   hmetspec<-draw(hsmet2)
  #   print(hmetspec)
  #   dev.off()
  #   
  # }
  

# 
# #
# # tSNE analysis
# #
# 
# if(tsne!=FALSE){
#   
# colnames(mat_hs_for_tsne3)<-gsub(make.names(colnames(mat_hs_for_tsne3),unique=T),pattern="\\.",replacement="-")
# 
# input_tsne<-rbind(mat_hs_for_tsne3,mat_meth3)
# input_tsne[is.na(input_tsne)]<-0
# 
# library(M3C)
# 
#        
# output_tsne<-paste(sapply(strsplit(output_hs,split=".png"),"[[",1),"tSNE.png",sep=".")
# 
# png(output_tsne,width=8,height=8)
# 
# print("create tsne")
# print(dim(input_tsne))
# 
# for(i in variables_to_use_tsne){
# 	
# 	print(i)
# 
# 	ann<-annotation2[,which(colnames(annotation2)%in%i)]
# 
# 	p<-tsne(input_tsne,labels=as.factor(ann), dotsize = 2,perplex=perplex,seed=123)
# 
# 	print(p)
# }
# 
# dev.off()
# }

}

