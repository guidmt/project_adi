library(GenomicRanges)
library(data.table)
library(tidyr)
library(methimpute)
library(snow)
library(BSgenome.Mmusculus.UCSC.mm10)
library(foreach)
library(annotatr)
library(parallel)
library(tidyverse)
library(GenomicRanges)

annots = c(
	   'mm10_genes_promoters',				               
	   'mm10_genes_exons')

annotations_order = build_annotations(genome = 'mm10', annotations = annots)
all_the_regions<-as.data.frame(annotations_order)[,c(1:3,6,10)]
all_the_regions$id<-sapply(strsplit(all_the_regions$id,split="\\:"),"[[",1)
all_the_regions<-unique(all_the_regions)

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis")
source("support_function_heatmap_mouse.ge.R")

# Load each file of methylation
input_met_mm_rid<-read.delim(file="/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/annotation_dmrs_adi_pval0.05_lfc1.txt",header=T)
input_met_mm_rid<-input_met_mm_rid[,c(1:3,19,17)]
colnames(input_met_mm_rid)<-c("chr","Cpg_start","Cpg_end","gene_symbol","status")
input_met_mm_rid[,1]<-gsub(input_met_mm_rid[,1],pattern="chr",replacement="")

df_ann<-unique(input_met_mm_rid)


setwd("/mnt/data/lab/gmt_data/data_brain/adi/CytosineReports")
list_meth_file<-dir()

#
# Define a function to get the methylation values
#

getSignalMeth<-function(X,meth2_chr,positions){
  
  start<-which.min(abs(positions - as.numeric(X[2])))
  end<-which.min(abs(positions - as.numeric(X[3])))
  
  meth2_chr2<-meth2_chr[start:end,]
  
  cov = meth2_chr2$count_methylated +  meth2_chr2$count_unmethylated
  idx_to_use<-which(cov>=5)
  meth2_chr2<-meth2_chr2[idx_to_use,]
  
  methSignal<- sum(meth2_chr2$count_methylated)/(sum(meth2_chr2$count_methylated)+sum(meth2_chr2$count_unmethylated))
  
  methSignal2<-paste(paste(unique(meth2_chr2[,1]),X[2],X[3],sep="_"),methSignal,sep=":")
  
}

#
# Store the results in some objects
#

list_samples<-vector(mode='list',length(list_meth_file))
df_meth_all<-data.frame()
signals_all_samples<-data.frame()

for(i in 1:length(list_meth_file)){
  
  current_name<-gsub(list_meth_file[i],pattern=".nsorted.deduplicated.CpG_report.txt.gz",replacement="")

  print(current_name)
  
  print("reading data")
  meth_table_current_mouse<-fread(file=list_meth_file[i])
  colnames(meth_table_current_mouse)<-c("chromosome","position","strand","count_methylated","count_unmethylated","c_context","trinucleotide_context")

  meth_table_current_mouse<-meth_table_current_mouse[,c(1,2,4,5)]
  setkey(meth_table_current_mouse,chromosome,position)
  
    list_chr<-grep(unique(all_the_regions[,1]),pattern="random",invert=T,value=T)
    
    all_chromosomes_signals<-data.frame()
    
    for(chr_region in list_chr){
    
    print(chr_region)
      
    meth2_chr<-meth_table_current_mouse[chromosome==chr_region]
    all_the_regions_chr<-all_the_regions[which(all_the_regions[,1]%in%chr_region),]

    n = 10
    split_df <- all_the_regions_chr %>% group_by(row_number() %/% n) %>% group_map(~ .x)
    
    #
    #Start parallel process
    #
    
    cl <- makeCluster(50, type = "SOCK") 
    
    clusterExport(cl, c("meth2_chr","getSignalMeth"), envir=environment())
    
    getStart<-parLapply(cl,split_df,fun=function(X){
      
      positions<-as.numeric(meth2_chr$position)
      
      input_list<-X
      
      meth2_chr<-meth2_chr
      
      res<-apply(input_list,1,getSignalMeth,meth2_chr=meth2_chr,positions)
  
    })
    
    stopCluster(cl)
    
    #
    # End parallel process
    #
    
    signals_for_regions<-cbind(samples=current_name,all_the_regions_chr,do.call(c,getStart))
  	
    all_chromosomes_signals<-rbind(all_chromosomes_signals,signals_for_regions)
    
    }
    
    signals_all_samples<-rbind(signals_all_samples,all_chromosomes_signals)
}

#
# Save the results
#

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/processing_meth_April")
save(signals_all_samples,file="Meth_processingGenomicFeatures.RData")

input_for_charts<-signals_all_samples[,c(1,5)]
input_for_charts$meth_signal<-as.numeric(sapply(strsplit(as.character(signals_all_samples[,7]),split="\\:"),"[[",2))
input_for_charts$status<-rep("NA",nrow(input_for_charts))
input_for_charts[grep(input_for_charts[,1],pattern=paste(c("wt","IWT"),collapse="|")),"status"]<-"WT"
input_for_charts[grep(input_for_charts[,1],pattern="IAP"),"status"]<-"MUT"
input_for_charts[grep(input_for_charts[,1],pattern="Control"),"status"]<-"CTRL"
input_for_charts$id<-factor(input_for_charts$id,level=c("promoter","exon"))

png("GlobalMethLevels.GenomicRegions.png", units="in", width=10, height=10, res=300)
p <- ggplot(input_for_charts, aes(x=samples, y=meth_signal,fill=status)) + 
  geom_boxplot()+scale_fill_manual(values = c("MUT" = "gold2", "WT" = "green4","CTRL"="dodgerblue2"))+ylab("M/(M+U)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~id)
print(p)
dev.off()

png("GlobalMethLevels.GenomicRegions.yesredo.nocontrol.png", units="in", width=10, height=10, res=300)
input_for_charts2<-input_for_charts[which(input_for_charts[,1]%in%c("IAP89","IAP89B","IAPA","IwtAP2R4a","IwtAP6c_redo","IWTAPB")),]
p <- ggplot(input_for_charts2, aes(x=samples, y=meth_signal,fill=status)) + 
  geom_boxplot()+scale_fill_manual(values = c("MUT" = "gold2", "WT" = "green4","CTRL"="dodgerblue2"))+ylab("M/(M+U)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~id)
print(p)
dev.off()

#
# Overlap with the DMR regions from the Canada People islands
#

input_met_mm<-read.delim(file="/mnt/data/lab/gmt_data/data_brain/adi/BriefSummaryReports/ListOfDifferentiallyMethylatedCGI-Mutants.txt",header=F)
input_met_mm_rid<-input_met_mm[,c(1:3)]
colnames(input_met_mm_rid)<-c("chr","start","end")

input_for_charts_cpg<-signals_all_samples[,c(2,3,4,1,5,7)]
input_for_charts_cpg$meth_signal<-as.numeric(sapply(strsplit(as.character(input_for_charts_cpg[,6]),split="\\:"),"[[",2))
input_for_charts_cpg$status<-rep("NA",nrow(input_for_charts_cpg))
input_for_charts_cpg[grep(input_for_charts_cpg$samples,pattern=paste(c("wt","IWT"),collapse="|")),"status"]<-"WT"
input_for_charts_cpg[grep(input_for_charts_cpg$samples,pattern="IAP"),"status"]<-"MUT"
input_for_charts_cpg[grep(input_for_charts_cpg$samples,pattern="Control"),"status"]<-"CTRL"
input_for_charts_cpg$id<-factor(input_for_charts_cpg$id,level=c("promoter","exon"))

input_met_mm_rid_gr<-makeGRangesFromDataFrame(input_met_mm_rid)
input_for_charts_cpg_gr<-makeGRangesFromDataFrame(input_for_charts_cpg,keep.extra.columns = TRUE)

cpg_signals_meth<-findOverlaps(input_met_mm_rid_gr,input_for_charts_cpg_gr,maxgap=1000)

qh<-queryHits(cpg_signals_meth)
sh<-subjectHits(cpg_signals_meth)

res_query<-data.frame(input_met_mm_rid[qh,])
res_subject<-data.frame(input_for_charts_cpg[sh,])

png("GlobalMethLevels.GenomicRegions.DMR_strict_Canada.png", units="in", width=10, height=10, res=300)
p <- ggplot(res_subject, aes(x=samples, y=meth_signal,fill=status)) + 
  geom_boxplot()+scale_fill_manual(values = c("MUT" = "gold2", "WT" = "green4","CTRL"="dodgerblue2"))+ylab("M/(M+U)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~id)
print(p)
dev.off()

png("GlobalMethLevels.GenomicRegions.DMR_strict_Canada.yesredo.nocontrol.png", units="in", width=10, height=10, res=300)
input_for_charts2<-res_subject[which(res_subject$samples%in%c("IAP89","IAP89B","IAPA","IwtAP2R4a","IwtAP6c_redo","IWTAPB")),]
p <- ggplot(input_for_charts2, aes(x=samples, y=meth_signal,fill=status)) + 
  geom_boxplot()+scale_fill_manual(values = c("MUT" = "gold2", "WT" = "green4","CTRL"="dodgerblue2"))+ylab("M/(M+U)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~id)
print(p)
dev.off()
