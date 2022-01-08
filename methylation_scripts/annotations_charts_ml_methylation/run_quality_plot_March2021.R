library(VennDetail)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggpubr)

source("plot_annotation_custom.R")

setwd("/mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methRegion/unfiltered/withLogFC")

unfiltered_with_logFC<-read.table(file="CpG_Islands_DMR-calls.txt",fill=TRUE,header=T)

setwd("/mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methRegion/unfiltered/DynmicFragmentation")

unfiltered_dynamicfragmentation<-read.table(file="CpG_Islands_DMR-calls.txt",fill=TRUE,header=T)

setwd("/mnt/data/lab/gmt_data/data_brain/meth_data_06_03_2020/adi/BriefSummaryReports")

DifferentialCGI<-read.table("ListOfDifferentiallyMethylatedCGI-Mutants.txt")
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
# Barchat with the % of regions that are methylated or hypo-methylated
#
matrix_DMRs<-unfiltered_with_logFC[,c(2:7)]


#
# Volcano plot mut vs wt 
#

pdf("volcano_plot_mut_vs_wt_regions_v1.pdf",width=5,height=5)

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


checkMethMatrix<-function(x){

	thr<-quantile(x)[2]
	print(thr)
        hyper<-length(x[x>thr])/length(x)
	hypo<-length(x[x<=thr])/length(x)
	res<-data.frame(hyper,hypo)

	return(res)
}

regions_df2<-unfiltered_with_logFC

mtx<-regions_df2[,c(2:7)]
mtx2<-do.call(rbind,apply(mtx,2,checkMethMatrix))

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")

pdf("DMRs_barplot_meth_perc_all.pdf")

barplot(t(mtx2), main="%Methylated regions",	
	xlab="Sample", col=c("red3","green3"),
	legend = rownames(t(mtx2)))

dev.off()

regions_df2<-unfiltered_with_logFC

status<-rep("no_sig",nrow(unfiltered_with_logFC))

status[which(regions_df2$pvalue<=0.05 & regions_df2$log2FoldChange>1)]<-"DMRs-up"
status[which(regions_df2$pvalue<=0.05 & regions_df2$log2FoldChange< -1)]<-"DMRs-down"

regions_df2$status<-status

mtx<-regions_df2[regions_df2$status!="no_sig",c(2:7)]
mtx2<-do.call(rbind,apply(mtx,2,checkMethMatrix))

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")

pdf("DMRs_barplot_meth_perc_p0.05_log2fc1.pdf")

barplot(t(mtx2), main="%Methylated regions",
	xlab="Sample", col=c("red3","green3"),
	legend = rownames(t(mtx2)))
dev.off()



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
				  quiet = FALSE,
				  minoverlap = 150)
#Why minoverlap 150? Is the best cut-off to reduce the number of regions but annotate each region that we get
#test<-data.frame(annotated_dmrs)
#dim(test)
#test2<-test[,c(17,27)]
#dim(unique(test2))
#test3<-unique(test2)
#table(test3[,1])
#DMRs-down   DMRs-up 
#      307      1744 

png("DMRs_sig_p0.05_lfc1.png",res=200,width=15,height=15,units="cm")
plot_annotation(
	    	annotated_regions = annotated_dmrs,
		plot_title = 'Dist. DMRs sig. p.value <= 0.05 and lfc 1')
dev.off()

pdf("DMRs_sig_p0.05_lfc1_piechart_up_down_December.pdf")
plot_annotationPieChart_GMT(annotated_regions = annotated_dmrs)
dev.off()

pdf("DMRs_sig_p0.05_lfc1_break_up_down_March2021_perc_selected.pdf")
plot_annotationGMT(annotated_regions = annotated_dmrs,perc=T,list_selection=c("promoters","introns","exons"),select=T)
dev.off()

pdf("DMRs_sig_p0.05_lfc1_break_up_down_March2021_perc_all.pdf")
plot_annotationGMT(annotated_regions = annotated_dmrs,perc=T,select=F)
dev.off()

#png("DMRs_sig_p0.05_lfc1_break_up_down_December.png",res=200,width=15,height=15,units="cm")
#plot_annotationGMT(annotated_regions = annotated_dmrs,perc=F)
#dev.off()

#
# Create the output for adi: a simplified version
#

res_annotation<-data.frame(annotated_dmrs)

#
# Create a Volcano Plot to color the regions
#


res_annotation_ffv<-res_annotation

# find which is the closest genomic regions associated with a CpG island
#

# Start CGI - end genomic annotation
res_annotation_ffv$start_meth_end_gr<-res_annotation$start-res_annotation$annot.end
res_annotation_ffv$start_gr_end_meth<-res_annotation$annot.start-res_annotation$end

res_annotation_ffv$mean_dist<-apply(res_annotation_ffv[,c(29,30)],1,mean)
res_annotation_ffv2<-unique(res_annotation_ffv)

res_annotation_ffv2$cpgIsland<-as.character(res_annotation_ffv2$cpgIsland)
list_cpg<-unique(res_annotation_ffv2$cpgIsland)

res_annotation_ffv3<-data.frame()

for(rsaf in 1:length(list_cpg)){
  
  col_to_use<-c("cpgIsland","start","end","annot.start","annot.end","mean_dist","status","annot.type","log2FoldChange","pval2")
  
  sub_raf2<-unique(res_annotation_ffv2[which(res_annotation_ffv2$cpgIsland%in%list_cpg[rsaf]),which(colnames(res_annotation_ffv2)%in%col_to_use)])
  
  if(nrow(sub_raf2[sub_raf2$annot.type%in%"mm10_genes_promoters",])==1){
    
    print(rsaf)
    
    if(abs(sub_raf2[sub_raf2$annot.type%in%"mm10_genes_promoters","mean_dist"])<=1000){
    
    print("TRUE") 
      
    res_annotation_ffv3<-rbind(res_annotation_ffv3,sub_raf2[sub_raf2$annot.type%in%"mm10_genes_promoters",])
    
    }
    
    }else{
    
  idx_raf_min<-which.max(sub_raf2$mean_dist)
  
  sub_raf3<-sub_raf2[idx_raf_min,]
  
  res_annotation_ffv3<-rbind(res_annotation_ffv3,sub_raf3)
  
  }
  
}



res_annotation_ffv3$new_annotation<-res_annotation_ffv3$annot.type
idx_regions_to_highlight<-which(res_annotation_ffv3$annot.type%in%c("mm10_genes_promoters","mm10_genes_exons"))

#
# This data.frame contains only the regions that are DMRs - This is a better representation of the reality respect the barchart
#

res_annotation_ffv3[-idx_regions_to_highlight,"new_annotation"]<-"other_genomics_DMRs"

input_volcano_rid_ann_DMRs<-res_annotation_ffv3[,c("cpgIsland","log2FoldChange","pval2","new_annotation","status")]
noSigRegions<-unique(regions_df[which(regions_df$status%in%c("no_sig")),c("cpgIsland","log2FoldChange","pval2","status")])
noSigRegions$new_annotation<-"other_genomics"

# Add regions that maybe have been lost during the annotation - 10 regions are lost during the annotation

input_volcano_tot_all_withann<-rbind(input_volcano_rid_ann_DMRs,noSigRegions)

idold<-unique(regions_df_filter$cpgIsland)
idnew<-unique(input_volcano_tot_all_withann$cpgIsland)
recovered_regions<-setdiff(idold,idnew)

recoveredDF<-unique(regions_df[which(regions_df$cpgIsland%in%recovered_regions),c("cpgIsland","log2FoldChange","pval2","status")])
recoveredDF$new_annotation<-"other_genomics"

input_volcano_tot_all_withann2<-rbind(input_volcano_tot_all_withann,recoveredDF)



pdf("VolcanoDMRs_colored_by_GenomicFeatures.pdf")

p<-ggplot(input_volcano_tot_all_withann2, aes(x=log2FoldChange, y=pval2)) +
  geom_point(aes(shape=new_annotation,color=new_annotation),size=0.6)+ 
  scale_color_manual(values=c("blue2","red2","gainsboro","lightblue4"))+labs(y= "-log(p.value)", x = "log2FoldChange")+ 
  theme_classic()+geom_vline(xintercept=-c(-1,0,1),linetype="longdash",color="gray2",size=0.2)+
  geom_hline(yintercept=-log(0.05,10),linetype="longdash",color="black",size=0.2)
print(p)

res<-table(input_volcano_tot_all_withann2$new_annotation,input_volcano_tot_all_withann2$status)

#This contigency table contains 1964 samples and not 1974 (delta = 10)
# DW 268 vs 271    (delta = 3)
# UP  1696 vs 1703 (delta = 7)

library(gridExtra)
print(grid.arrange(tableGrob(res)))

dev.off()


res_annotation2<-unique(res_annotation[,c(1:18,27,ncol(res_annotation))])
write.table(res_annotation2,file="annotation_dmrs_adi_pval0.05_lfc1.txt",sep="\t",row.names=F,quote=F)

save_clean_annotation<-merge(input_volcano_tot_all_withann2,unique(res_annotation2[,c("cpgIsland","annot.symbol")]),by="cpgIsland",all.x=T)
checkDR<-data.frame(table(save_clean_annotation[,1]))
ID<-unique(checkDR[checkDR[,2]>=2,1])
idx_clean<-which(save_clean_annotation[,1]%in%ID & is.na(save_clean_annotation$annot.symbol))
save_clean_annotation2<-save_clean_annotation[-idx_clean,]

write.table(save_clean_annotation2,file="regions_dmrs_adi_pval0.05_lfc1.GenomicRegions_and_Genes.Aug.txt",sep="\t",row.names=F,quote=F)

#
# add the information about the TFs related with the genes
#
library(readxl)

read_go_tf_Ray_up<-read_excel("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/DMR_enrichR.v3.TF.byType.xlsx","DMRs-up")

read_go_tf_Ray_dw<-read_excel("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/DMR_enrichR.v3.TF.byType.xlsx","DMRs-down")
all_go_Ray<-data.frame(rbind(read_go_tf_Ray_up,read_go_tf_Ray_dw))
#I do in this way, because i want a unique association between TFs and Genes - indipendently by the species (MM/HS)
all_go_Ray$Term<-sapply(strsplit(all_go_Ray$Term,split=" "),"[[",1)

  
matrixTFforGenes<-lapply(split(all_go_Ray,all_go_Ray$Term),FUN=function(X){
  
  unique_genes<-unique(unlist(strsplit(X$Genes,split="\\;")))
  unique_TFs<-unique(X$Term)
  
  res<-data.frame(TF=rep(unique_TFs,length(unique_genes)),genes=unique_genes)
  
})
matrixTFforGenes2<-unique(do.call(rbind,matrixTFforGenes))
genes_for_tf_df<-as.data.frame.matrix(table(matrixTFforGenes2$genes,matrixTFforGenes2$TF))
genes_for_tf_df2<-data.frame(Genes=rownames(genes_for_tf_df),genes_for_tf_df)

require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
humanx = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes_for_tf_df2$Genes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

#there is not  PRC and TBX2 (BUT TBX20 that is selected in the grep inside the next loop for)
list_TFs<-c("Genes","TP53","KDM2B","OLIG2","KDM4A","SUZ12","EZH2","DMRT1","GATA1","GATA2","GATA3","GATA4","EOMES","FOXA1","RUNX3","SPIB","TFAP2A","MYC","MYCN","E2F1","E2F2","E2F4")

genes_for_tf_df3<-merge(humanx,genes_for_tf_df2[,which(colnames(genes_for_tf_df2)%in%list_TFs)],by.x="HGNC.symbol",by.y="Genes")

#
# Merge the annotation table with the TFs
#
save_clean_annotation3<-merge(save_clean_annotation2,genes_for_tf_df3,by.x="annot.symbol",by.y="MGI.symbol",all.x=T,all.y=T)
save_clean_annotation3[is.na(save_clean_annotation3)]<-0

write.table(save_clean_annotation3,file="regions_dmrs_adi_pval0.05_lfc1.GenomicRegions_and_Genes.Aug.WithTFs.txt",sep="\t",row.names=F,quote=F)


regions_to_save<-unique(res_annotation2[,c(1,2,3)])
write.table(regions_to_save,file="regions_dmrs_adi_pval0.05_lfc1.REGIONS.bed",sep="\t",row.names=F,quote=F)

regions_to_save_pe<-unique(res_annotation2[which(res_annotation2$annot.type%in%c("mm10_genes_promoters","mm10_genes_exons")),c(1,2,3)])
write.table(regions_to_save_pe,file="regions_dmrs_adi_pval0.05_lfc1.PromotersExons.bed",sep="\t",row.names=F,quote=F)

write.table(unique(cbind(res_annotation2[,c(1,2)],res_annotation2[,c(2)]+1)),file="regions_dmrs_adi_pval0.05_lfc1.START.bed",col.names=F,sep="\t",row.names=F,quote=F)


# UP
write.table(unique(res_annotation2[res_annotation2$status=="DMRs-up",c(1:3)]),file="regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.bed",col.names=F,sep="\t",row.names=F,quote=F)
upDF<-unique(res_annotation2[which(res_annotation2$status%in%"DMRs-up" & res_annotation2$annot.type %in% c("mm10_genes_promoters","mm10_genes_exons")),c(1:3)])
write.table(upDF,file="regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.bed",col.names=F,sep="\t",row.names=F,quote=F)

# DOWN
write.table(unique(res_annotation2[res_annotation2$status=="DMRs-down",c(1:3)]),file="regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.bed",col.names=F,sep="\t",row.names=F,quote=F)
dwDF<-unique(res_annotation2[which(res_annotation2$status%in%"DMRs-down" & res_annotation2$annot.type %in% c("mm10_genes_promoters","mm10_genes_exons")) ,c(1:3)])
write.table(dwDF,file="regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.bed",col.names=F,sep="\t",row.names=F,quote=F)


# UP Promoters
# write.table(unique(res_annotation2[res_annotation2$status=="DMRs-up" & res_annotation2$annot.type=="mm10_genes_promoters",c(1:3)]),file="regions_dmrs_adi_pval0.05_lfc1.REGIONS.Promoters.UP.bed",col.names=F,sep="\t",row.names=F,quote=F)

# DOWN Prometers
# write.table(unique(res_annotation2[res_annotation2$status=="DMRs-down" & res_annotation2$annot.type=="mm10_genes_promoters",c(1:3)]),file="regions_dmrs_adi_pval0.05_lfc1.REGIONS.Promoters.DOWN.bed",col.names=F,sep="\t",row.names=F,quote=F)

library(readxl)

read_go_tf_Ray_up<-read_excel("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/DMR_enrichR.v3.TF.byType.xlsx","DMRs-up")

read_go_tf_Ray_dw<-read_excel("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/DMR_enrichR.v3.TF.byType.xlsx","DMRs-down")
all_go_Ray<-data.frame(rbind(read_go_tf_Ray_up,read_go_tf_Ray_dw))

#idx_mouse<-grep(all_go_Ray$Term,pattern="mouse")
#idx_mouse2<-grep(all_go_Ray$Term,pattern="Mouse")
#idx_mm9<-grep(all_go_Ray$Term,pattern="mm9")


#all_go_Ray2<-all_go_Ray[c(idx_mouse,idx_mouse2,idx_mm9),]
all_go_Ray2<-all_go_Ray
#Keep space after MYC for grep and distinguish by MYCN
#DMTR
#there is not  PRC and TBX2 (BUT TBX20 that is selected in the grep inside the next loop for)
list_TFs<-paste("^",c("TP53","PRC2","KDM2B","OLIG2","KDM4A","SUZ12","EZH2","DMRT1","GATA1","GATA2","GATA3","GATA4","EOMES","FOXA1","RUNX3","SPIB","TBX2","TFAP2A","MYC ","MYCN","E2F1","E2F2","E2F4"),sep="")

regions_to_save2<-res_annotation2[,c(1,2,3,19,20)]

require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
humanx = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = regions_to_save2$annot.symbol , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

regions_for_tf<-merge(regions_to_save2,humanx,by.x="annot.symbol",by.y="MGI.symbol")

for(i in list_TFs){

	tf_up_and_down<- unique(unlist(strsplit(all_go_Ray2[grep(all_go_Ray2$Term,pattern=i),10],split="\\;")))
	regions_for_tf2<-unique(regions_for_tf[which(regions_for_tf$HGNC.symbol%in%tf_up_and_down),c(2:6)])
	regions_for_tf2_pe<-unique(regions_for_tf[which(regions_for_tf$HGNC.symbol%in%tf_up_and_down & regions_for_tf$annot.type %in% c("mm10_genes_promoters","mm10_genes_exons")),c(2:6)])

	i2<-gsub(gsub(i,pattern="\\^",replacement=""),pattern=" ",replacement="")
	output_file<-paste(i2,"regions_dmrs_adi_pval0.05_lfc1.REGIONS.txt",sep="_")
        output_file2<-paste(i2,"regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt",sep="_")
	write.table(regions_for_tf2[,-c(4:5)],file=output_file,sep="\t",row.names=F,quote=F,col.names=F)
        write.table(regions_for_tf2[,-c(4:5)],file=output_file2,sep="\t",row.names=F,quote=F,col.names=F)

	all_go_Ray2_up<-all_go_Ray2[all_go_Ray2[,1]%in%"DMRs-up",]
	tf_up<- unique(unlist(strsplit(all_go_Ray2_up[grep(all_go_Ray2_up$Term,pattern=i),10],split="\\;")))	
	regions_for_tf2_up<-unique(regions_for_tf[which(regions_for_tf$HGNC.symbol%in%tf_up),c(2:4)])
	regions_for_tf2_up_pe<-unique(regions_for_tf[which(regions_for_tf$HGNC.symbol%in%tf_up & regions_for_tf$annot.type %in% c("mm10_genes_promoters","mm10_genes_exons")),c(2:4)])

	i2<-gsub(gsub(i,pattern="\\^",replacement=""),pattern=" ",replacement="")
        output_file<-paste(i2,"regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.txt",sep="_")	
	output_file2<-paste(i2,"regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt",sep="_")
	write.table(regions_for_tf2_up,file=output_file,sep="\t",row.names=F,quote=F,col.names=F)
	write.table(regions_for_tf2_up_pe,file=output_file2,sep="\t",row.names=F,quote=F,col.names=F)

	all_go_Ray2_dw<-all_go_Ray2[all_go_Ray2[,1]%in%"DMRs-down",]
	tf_dw<- unique(unlist(strsplit(all_go_Ray2_dw[grep(all_go_Ray2_dw$Term,pattern=i),10],split="\\;")))   
	regions_for_tf2_dw<-unique(regions_for_tf[which(regions_for_tf$HGNC.symbol%in%tf_dw),c(2:4)])
        regions_for_tf2_dw_pe<-unique(regions_for_tf[which(regions_for_tf$HGNC.symbol%in%tf_dw & regions_for_tf$annot.type %in% c("mm10_genes_promoters","mm10_genes_exons")),c(2:4)])

	i2<-gsub(gsub(i,pattern="\\^",replacement=""),pattern=" ",replacement="")
	output_file<-paste(i2,"regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.txt",sep="_")
	output_file2<-paste(i2,"regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt",sep="_")
	write.table(regions_for_tf2_dw,file=output_file,sep="\t",row.names=F,quote=F,col.names=F)
        write.table(regions_for_tf2_dw_pe,file=output_file2,sep="\t",row.names=F,quote=F,col.names=F)

}



#
# Overlap with gene-list C1A vs C1B
#

skip=TRUE

if(skip!=TRUE){

a<-c("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_17_08_2020/res_17_08_2020/DESeq2_rmv5_CSHC1A_Run1C2_Adi18//Comp2_C1A_vs_C1B_0.05_2.ort.txt")

list_file_c1a_c1b<-c(a)
names(list_file_c1a_c1b)<-c("DESeq2_rmv5_Adi18_10_22_CSHS15")

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

}

#
# Overlap with gene-list C1 vs C2
#

a<-c("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_11_11_2020/res_11_11_2020/DESeq2_rmv5_Adi18_10_22_CSHS15/Comp1_C1_vs_C2/Comp1_C1_vs_C2_0.01_2.ort.txt")

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")

list_file_c1a_c1b<-c(a)
names(list_file_c1a_c1b)<-c("DESeq2_rmv5_Adi18_10_22_CSHS15")

for(i in 1:length(list_file_c1a_c1b)){

	tab_input<-read.delim(file=list_file_c1a_c1b[i],stringsAsFactors=F)
	out2<-gsub(sapply(strsplit(list_file_c1a_c1b[i],split="/"),"[[",13),pattern=".txt",replacement="")

        string_output<-names(list_file_c1a_c1b)[i]

        gene_lists<-split(tab_input,tab_input$status)

  	gene_lists_down<-gene_lists[["down"]][,2]

	gene_lists_up<-gene_lists[["up"]][,2]


	library(VennDetail)
	
        dmrs_up=unique(res_annotation[res_annotation$status%in%"DMRs-up","annot.symbol"])
	dmrs_down=unique(res_annotation[res_annotation$status%in%"DMRs-down","annot.symbol"])

	dmrs_up<-dmrs_up[!is.na(dmrs_up)]	
	dmrs_down<-dmrs_down[!is.na(dmrs_down)]


	ven <- venndetail(list(DMRs_up = dmrs_up,
			       DMRs_down = dmrs_down,
			       down = gene_lists_down,			
			       up = gene_lists_up))


	venn_overlap<-ven@wide
	write.table(venn_overlap,paste(string_output,".",out2,".ovp_dmr_pvalue0.05.lfc1_December.txt",sep=""),sep="\t")

	pdf(paste(string_output,".",out2,".ovp_dmr_pvalue0.05.lfc1_December.pdf",sep=""))	
	p<-plot(ven,type = "upset")
	print(p)
	dev.off()

}


#
# Integrate gene-expression with methylation data
#

source("function_integrate_meth_ge.R")

methWithGE(methylation_data=res_annotation2,current_qval=0.01,current_lfc=2,string_out="December_meth_pvalue0.05_lfc1")
