methWithGE<-function(methylation_data,current_qval,current_lfc,lfc_meth=1,string_out){
	
		methylation_data_uni2<-unique(methylation_data[,c(1,2,3,13,19)])
		
		list_dir<-c("DESeq2_rmv5_Adi18_10_22_CSHS15")
		
		basic_dir<-c("/mnt/data/lab/gmt_data/data_brain/adi_data_november")

		for(mydir in 1:length(list_dir)){
			
			mydir2<-list_dir[mydir]
			
			print(mydir2)

		        input_file<-paste(paste(basic_dir,mydir2,"Comp1_C1_vs_C2",sep="/"),"IDH.Comp1_C1_vs_C2.DESeq2.sigDE.csv",sep="/")

		        degs_in_mouse<-read.csv(file=input_file,stringsAsFactors=F)

		        statistics_degs<-degs_in_mouse[,c(1:8)]
			
		      	selected_genes_with_degs<-statistics_degs[statistics_degs$padj<=current_qval & abs(statistics_degs$log2FoldChange)>=current_lfc,]

			      common_genes_meth_ge<-intersect(selected_genes_with_degs[,1],methylation_data_uni2[,5])

      			input_ge<-statistics_degs[statistics_degs[,1]%in%common_genes_meth_ge,c(1,4)]
      			colnames(input_ge)[1:2]<-c("genes_symbol","log2fc_ge")
      			
      			input_meth<-unique(methylation_data_uni2[methylation_data_uni2$annot.symbol%in%common_genes_meth_ge,c(5,4)])
      			colnames(input_meth)[1:2]<-c("genes_symbol","log2fc_meth")
      		
      			input_meth2<-aggregate(input_meth[,2], by = list(input_meth$genes_symbol), mean)
      			colnames(input_meth2)[1:2]<-c("genes_symbol","log2fc_meth")
      
      			all_data<-merge(input_ge,input_meth2,by="genes_symbol")
      
      			all_data$status<-rep("no_sig",nrow(all_data))
      
      			all_data[all_data$log2fc_ge >0 & all_data$log2fc_meth >0,"status"] <- "up-ge/up-meth"			
      			
      			all_data[all_data$log2fc_ge <0 & all_data$log2fc_meth >0,"status"] <- "down-ge/up-meth"
      	 
       			all_data[all_data$log2fc_ge >0 & all_data$log2fc_meth <0,"status"] <- "up-ge/down-meth"
      		
      			all_data[all_data$log2fc_ge <0 & all_data$log2fc_meth <0,"status"] <- "down-ge/down-meth"
      
      			write.table(all_data,paste("INT_GE_ME",paste(mydir2,current_qval,current_lfc,sep="_"),string_out,"txt",sep="."),sep="\t",row.names=F,quote=F)
      			pdf(paste("INT_GE_ME",paste(mydir2,current_qval,current_lfc,sep="_"),string_out,"pdf",sep="."))
            
      			require(gg.gap)
      			
      			p<-ggplot(all_data, aes(x=log2fc_meth, y=log2fc_ge)) +
      			   geom_point(aes(color=status, fill=status,alpha=0.5),size=3)+scale_shape_manual(values=c(20, 20,1))+scale_color_manual(values=c("blue2","red2","green2", "orange2"))+labs(y="log2FoldChange (ge)", x = "log2FoldChange (meth)") + theme_classic()+geom_vline(xintercept=-c(-lfc_meth,lfc_meth),linetype="longdash",color="black",size=0.2)+geom_hline(yintercept=c(-current_lfc,current_lfc),linetype="longdash",color="black",size=0.2)
            
      			print(p)
            
      			# #https://github.com/ChrisLou-bioinfo/gg.gap
      			# 
      			# print(p %>% gg.gap(ylim=c(min(all_data$log2fc_ge)-1,max(all_data$log2fc_ge))+1,segments = list(c(-1,1)),tick_width = c(0.5,0.5)))
      			require(plotrix)
      			gap.plot(all_data$log2fc_meth,
      			         all_data$log2fc_ge,
      			         gap=c(-1.5,1.5), 
      			         xlab="log2FoldChange (meth)",
      			         ylab="log2FoldChange (ge)",
      			         gap.axis="y",
      			         ytics=c(-6,-4,-2,-1,0,1,2,4,6),
      			         xtics=c(-6,-4,-2,-1,0,1,2,4,6),
      			         pch=19,
      			         breakcol = "white",
      			         col= as.numeric(as.factor(all_data$status)))
      			abline(h=c(-2,-1),lty =2)
      			# 2 is the axis on the left, -1.5 the position of breaks
      			axis.break(2,-1.5,style="slash")
      			axis.break(4,-1.5,style="slash")
      			# abline(h=c(-2,1),lty =2)
      			abline(v=c(-lfc_meth,lfc_meth),lty =2)
      			
      			dev.off()

		}
}
