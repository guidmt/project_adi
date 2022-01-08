convertMouseGeneList <- function(x){
		          
		  require("biomaRt")
		  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
		  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
		          
		  humanx = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

		  return(humanx)
			            
}


convertMouseGeneList2 <- function(x){

	require("biomaRt")
	human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",version="102")
	mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",version="102")
		                            
	humanx = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

	return(humanx)
				                                        
}


integrated_conversion_mouse_human<-function(mgi,ensembl){

               #common elements between mouse and human 
               mgi_table<-mgi[,c(1,5)]
                                               
               specifics_mgi<-setdiff(x=mgi_table[,2],y=ensembl[,1])
               specifics_mgi_table<-mgi_table[mgi_table[,2]%in%specifics_mgi,]
               colnames(specifics_mgi_table)<-c("hs","mm")

               specifics_ensembl<-setdiff(x=ensembl[,1],y=mgi_table[,2])
               specifics_ensembl_table<-ensembl[ensembl[,1]%in%specifics_ensembl,c(2,1)]
               colnames(specifics_ensembl_table)<-c("hs","mm")

               #take the common genes from mgi
               common_mgi_ensembl<-intersect(mgi_table[,2],ensembl[,1])

               common_table<-mgi_table[mgi_table[,2]%in%common_mgi_ensembl,]
               colnames(common_table)<-c("hs","mm")

               res_annotation<-rbind(common_table,specifics_ensembl_table,specifics_mgi_table)
	
	       return(res_annotation)
}
