library(VennDetail)
library(ggplot2)
source("plot_annotation_custom.R")

setwd("/mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methRegion/unfiltered/withLogFC")

unfiltered_with_logFC<-read.table(file="CpG_Islands_DMR-calls.txt",fill=TRUE,header=T)

setwd("/mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methRegion/unfiltered/DynmicFragmentation")

unfiltered_dynamicfragmentation<-read.table(file="CpG_Islands_DMR-calls.txt",fill=TRUE,header=T)

setwd("/mnt/data/lab/gmt_data/data_brain/meth_data_06_03_2020/adi/BriefSummaryReports")

DifferentialCGI<-read.table("ListOfDifferentiallyMethylatedCGI-Mutants.txt")
write.table(DifferentialCGI[,c(1:3)],file="ListOfDifferentiallyMethylatedCGI-Mutants.bed",sep="\t",quote=F,row.names=F,col.names=F)

DifferentialCGI$ID<-paste(DifferentialCGI[,1],DifferentialCGI[,2],DifferentialCGI[,3],sep="_")

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")

#
# Venn diagram to check
#

ven <- venndetail(list(unfiltered_with_logFC = unfiltered_with_logFC [,1], 
		       unfiltered_dynamicfragmentation = unfiltered_dynamicfragmentation[,1],		      
   		       DifferentialCGI = DifferentialCGI$ID))

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")

png("Adi_comparison_dmr_regions_ven.png")
plot(ven)
dev.off()

png("Adi_comparison_dmr_regions_upset.png")
plot(ven,type="upset")
dev.off()

#
# Volcano plot mut vs wt 
#

png("volcano_plot_mut_vs_wt_regions_v1.png",res=200,width=15,height=15,units="cm")

regions_df<-unfiltered_with_logFC

status<-rep("no_sig",nrow(unfiltered_with_logFC))

status[which(regions_df$pvalue<=0.05 & regions_df$log2FoldChange>1)]<-"DMRs-up"
status[which(regions_df$pvalue<=0.05 & regions_df$log2FoldChange< -1)]<-"DMRs-down"

regions_df$status<-status
regions_df$pval2<-(-log(regions_df$pvalue,10))

p<-ggplot(regions_df, aes(x=log2FoldChange, y=pval2)) +
  geom_point(aes(shape=status,color=status),size=0.6)+scale_shape_manual(values=c(20, 20, 1))+ scale_color_manual(values=c("blue2","red2", "gainsboro"))+labs(y= "-log(p.value)", x = "log2FoldChange")+ theme_classic()+geom_vline(xintercept=-c(-1,0,1),linetype="longdash",color="black",size=0.2)+geom_hline(yintercept=-log(0.05,10),linetype="longdash",color="black",size=0.2)

print(p)

dev.off()

table(regions_df$status,regions_df$DynamicallyDetected)

nrow(regions_df[regions_df$DynamicallyDetected=="yes" & regions_df$status=="DMRs-up",])


#
# Start the annotation and the integration with the gene-expression values C1A vs C1B
#

regions_df_filter<-regions_df[regions_df$status!="no_sig",]
significant_regions<-do.call(rbind,strsplit(as.character(regions_df_filter[,1]),split="\\_"))

start<-as.integer(significant_regions[,2])
end<-as.integer(significant_regions[,3])

significant_regions2<-data.frame(chr=significant_regions[,1],start=start,end=end)

regions_df_filter2<-cbind(significant_regions2,regions_df_filter)

library(annotatr)

annots = c(
	   'mm10_genes_1to5kb',
	   'mm10_genes_promoters',
	   'mm10_genes_5UTRs',
	   'mm10_genes_exons',
	   'mm10_genes_intronexonboundaries',
	   'mm10_genes_introns',
	   'mm10_genes_3UTRs',
	   'mm10_genes_intergenic')

annotations_order = build_annotations(genome = 'mm10', annotations = annots)

significant_regions2_gr<-makeGRangesFromDataFrame(regions_df_filter2,keep.extra.columns=T)

annotated_dmrs = annotate_regions(
				  regions = significant_regions2_gr,
				  annotations = annotations_order,
				  ignore.strand = TRUE,
				  quiet = FALSE)

png("DMRs_sig_p0.05_lfc1.png",res=200,width=15,height=15,units="cm")
plot_annotation(
	    	annotated_regions = annotated_dmrs,
		plot_title = 'Dist. DMRs sig. p.value <= 0.05 and lfc 1')
dev.off()

png("DMRs_sig_p0.05_lfc1_break_up_down.png",res=200,width=15,height=15,units="cm")
plot_annotationGMT(annotated_regions = annotated_dmrs)
dev.off()

#
# Create the output for adi: a simplified version
#

res_annotation<-data.frame(annotated_dmrs)
res_annotation2<-unique(res_annotation[,c(1:18,27,ncol(res_annotation))])

write.table(res_annotation2,file="annotation_dmrs_adi_pval0.05_lfc1.txt",sep="\t",row.names=F,quote=F)

regions_to_save<-unique(res_annotation2[,c(1,2,3,19)])

write.table(regions_to_save,file="regions_dmrs_adi_pval0.05_lfc1.txt",sep="\t",row.names=F,quote=F)


#
# Overlap with gene-list C1A vs C1B
#


a<-c("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_17_08_2020/res_17_08_2020/DESeq2_rmv5_CSHC1A_Run1C2_Adi18/Comp2_C1A_vs_C1B/Comp2_C1A_vs_C1B_0.05_2.ort.txt")

b<-c("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_17_08_2020/res_17_08_2020/DESeq2_rmv5_rmvAdi18/Comp2_C1A_vs_C1B/Comp2_C1A_vs_C1B_0.05_2.ort.txt")

c1<-c("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_17_08_2020/res_17_08_2020/DESeq2_rmv5_rmvCSHC1A_rmvAdi18/Comp2_C1A_vs_C1B/Comp2_C1A_vs_C1B_0.05_2.ort.txt")

list_file_c1a_c1b<-c(a,b,c1)
names(list_file_c1a_c1b)<-c("DESeq2_rmv5_CSHC1A_Run1C2_Adi18","DESeq2_rmv5_rmvAdi18","DESeq2_rmv5_rmvCSHC1A_rmvAdi18")

for(i in 1:length(list_file_c1a_c1b)){

	tab_input<-read.delim(file=list_file_c1a_c1b[i],stringsAsFactors=F)
	out2<-gsub(sapply(strsplit(list_file_c1a_c1b[i],split="/"),"[[",13),pattern=".txt",replacement="")

	string_output<-names(list_file_c1a_c1b)[i]

	gene_lists<-split(tab_input,tab_input$status)

	gene_lists_down<-gene_lists[["down"]][,2]
        gene_lists_up<-gene_lists[["up"]][,2]

	library(VennDetail)

	dmrs=unique(res_annotation[res_annotation$status%in%"DMRs-up","annot.symbol"])

	dmrs<-dmrs[!is.na(dmrs)]

	ven <- venndetail(list(DMRs= dmrs,
			       down = gene_lists_down,
			       up = gene_lists_up))

        venn_overlap<-ven@wide

	write.table(venn_overlap,paste(string_output,".",out2,".ovp_dmr_pvalue0.05.lfc1.txt",sep=""),sep="\t")

	pdf(paste(string_output,".",out2,".ovp_dmr_pvalue0.05.lfc1.pdf",sep=""))
	p<-plot(ven,type = "upset")
	print(p)
	dev.off()

}

#
# Overlap with gene-list C1 vs C2
#

a<-c("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_17_08_2020/res_17_08_2020/DESeq2_rmv5_CSHC1A_Run1C2_Adi18/Comp1_C1_vs_C2/Comp1_C1_vs_C2_0.001_2.ort.txt")

b<-c("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_17_08_2020/res_17_08_2020/DESeq2_rmv5_rmvAdi18/Comp1_C1_vs_C2/Comp1_C1_vs_C2_0.001_2.ort.txt")

c1<-c("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_17_08_2020/res_17_08_2020/DESeq2_rmv5_rmvCSHC1A_rmvAdi18/Comp1_C1_vs_C2/Comp1_C1_vs_C2_0.001_2.ort.txt")

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")

list_file_c1a_c1b<-c(a,b,c1)
names(list_file_c1a_c1b)<-c("DESeq2_rmv5_CSHC1A_Run1C2_Adi18","DESeq2_rmv5_rmvAdi18","DESeq2_rmv5_rmvCSHC1A_rmvAdi18")

for(i in 1:length(list_file_c1a_c1b)){

	tab_input<-read.delim(file=list_file_c1a_c1b[i],stringsAsFactors=F)
	out2<-gsub(sapply(strsplit(list_file_c1a_c1b[i],split="/"),"[[",13),pattern=".txt",replacement="")

        string_output<-names(list_file_c1a_c1b)[i]

        gene_lists<-split(tab_input,tab_input$status)

  	gene_lists_down<-gene_lists[["down"]][,2]

	gene_lists_up<-gene_lists[["up"]][,2]


	library(VennDetail)
	
        dmrs=unique(res_annotation[res_annotation$status%in%"DMRs-up","annot.symbol"])

	dmrs<-dmrs[!is.na(dmrs)]

	ven <- venndetail(list(DMRs= dmrs,
			       down = gene_lists_down,			
			       up = gene_lists_up))


	venn_overlap<-ven@wide
	write.table(venn_overlap,paste(string_output,".",out2,".ovp_dmr_pvalue0.05.lfc1.txt",sep=""),sep="\t")

	pdf(paste(string_output,".",out2,".ovp_dmr_pvalue0.05.lfc1.pdf",sep=""))	
	p<-plot(ven,type = "upset")
	print(p)
	dev.off()

}


#
# Integrate gene-expression with methylation data
#

source("function_integrate_meth_ge.R")

methWithGE(methylation_data=res_annotation2,current_qval=0.001,current_lfc=2,string_out="meth_pvalue0.05_lfc1")
