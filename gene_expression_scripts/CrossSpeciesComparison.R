CrossSpeciesComparison<-function(mat_exp,ortholog_table_res,species="hs"){
	  
    zscore=function(x){(x-mean(x))/sd(x)}
  
  	mat_exp_temp<-log(mat_exp+1,2)

    # Get the common genes between mm and hs
    if(species=="hs"){
	      mat_for_score<-t(apply(mat_exp_temp[rownames(mat_exp_temp)%in%ortholog_table_res[,1],],1,zscore))
        }else{
	      mat_for_score<-t(apply(mat_exp_temp[rownames(mat_exp_temp)%in%ortholog_table_res[,2],],1,zscore)) 
      
	      }
    
    if(species=="hs"){
	         
	  up_genes<-ortholog_table_res[ortholog_table_res$status=="up",1]
	  down_genes<-ortholog_table_res[ortholog_table_res$status=="down",1]
            
    }else{
		 
	  up_genes<-ortholog_table_res[ortholog_table_res$status=="up",2]
    down_genes<-ortholog_table_res[ortholog_table_res$status=="down",2]
	  
     }
      
      scores_up<-apply(mat_for_score[rownames(mat_for_score)%in%up_genes,],2,mean,na.rm=T)
      scores_down<-apply(mat_for_score[rownames(mat_for_score)%in%down_genes,],2,mean,na.rm=T)
      
      scores_global<-scores_up-scores_down
      
      # Create a second score based on gsva  
      gsva_up <- as.numeric(gsva(as.matrix(mat_exp_temp), list(up_genes),min.sz=5, verbose=TRUE))
      gsva_down<- as.numeric(gsva(as.matrix(mat_exp_temp), list(down_genes),min.sz=5, verbose=TRUE))

      gsva_scores<-gsva_up-gsva_down
      
      df_scores<-data.frame(ID=names(scores_global),scores_global,gsva_scores)
	  	  
      return(df_scores)
	  
}
