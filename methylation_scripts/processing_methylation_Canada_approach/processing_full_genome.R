library(GenomicRanges)
library(data.table)
library(tidyr)
library(methimpute)
library(snow)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(reshape)

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

#Read each file
list_samples<-vector(mode='list',length(list_meth_file))
df_meth_all<-data.frame()

for(i in 1:length(list_meth_file)){
  
  print(i)
  
  current_name<-gsub(list_meth_file[i],pattern=".nsorted.deduplicated.CpG_report.txt.gz",replacement="")
    
  print("reading data")
  meth_table_current_mouse<-data.frame(fread(file=list_meth_file[i]))
  colnames(meth_table_current_mouse)<-c("chromosome","position","strand","count_methylated","count_unmethylated","c_context","trinucleotide_context")
  
  meth2<-meth_table_current_mouse
  cov = meth2$count_methylated +  meth2$count_unmethylated
  idx_to_use<-which(cov>=5)
  meth2<-meth2[idx_to_use,]
  
  list_meth_chr<-split(meth2,meth2$chromosome)
  
  #define a function to compute the methylation values
  computeMeth<-function(Xbin){
    
    # cov  <- Xbin$count_methylated +  Xbin$count_unmethylated
    # Xcov <- cbind(Xbin,cov)
    
    methSignal<- (Xbin$count_methylated)/(Xbin$count_methylated+Xbin$count_unmethylated)
    
    methSignal2<- cbind(Xbin,methSignal)
    
    return(methSignal2)
    
  }
    
  #One position of the genome correspond with 1bp
  split_genome<-function(X,binsize=1000){
  
  df.new = split(X, (seq(nrow(X))-1) %/% binsize) 
  
  res<-do.call(rbind,lapply(df.new,computeMeth))
  
  }
  
  cl <- makeCluster(60) 
  clusterExport(cl, "computeMeth")

  rsp_temp<-parLapply(cl, list_meth_chr, split_genome)

  resSignalProcessing<-do.call(rbind,rsp_temp)
  
  stopCluster(cl)

  df_meth<-data.frame(current_name,as.numeric(resSignalProcessing$methSignal))
  df_meth_all<-rbind(df_meth_all,df_meth)
  
}

colnames(df_meth_all)[c(1:2)]<-c("samples","meth_signals")

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/processing_meth_April")
save(df_meth_all,file="Meth_processingJuly.1000kb.RData")


df_meth_all$status<-rep("NA",nrow(df_meth_all))
df_meth_all[grep(df_meth_all[,1],pattern=paste(c("wt","IWT"),collapse="|")),"status"]<-"WT"
df_meth_all[grep(df_meth_all[,1],pattern="IAP"),"status"]<-"MUT"
df_meth_all[grep(df_meth_all[,1],pattern="Control"),"status"]<-"CTRL"

png("GlobalMethLevels.1000kb.png", units="in", width=10, height=10, res=300)
p <- ggplot(df_meth_all, aes(x=samples, y=meth_signals,fill=status)) + 
	geom_boxplot()+scale_fill_manual(values = c("MUT" = "gold2", "WT" = "green4","CTRL"="dodgerblue2"))+ylab("M/(M+U)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
dev.off()

png("GlobalMethLevels.1000kb.yesredo.nocontrol.png", units="in", width=10, height=10, res=300)
df_meth_all2<-df_meth_all[df_meth_all[,1]%in%c("IAP89","IAP89B","IAPA","IwtAP2R4a","IwtAP6c_redo","IWTAPB"),]
p <- ggplot(df_meth_all2, aes(x=samples, y=meth_signals,fill=status)) + 
  geom_boxplot()+scale_fill_manual(values = c("MUT" = "gold2", "WT" = "green4","CTRL"="dodgerblue2"))+ylab("M/(M+U)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
dev.off()

#
# Barplot methylation at sample levels
#

cut_off_methylation<-0.7
  
list_samples<-unique(df_meth_all2[,1])

stats_meth<-data.frame()

for(lss in list_samples){
  
  subset_df<-df_meth_all2[which(df_meth_all2[,1]%in%lss),]
  status_meth<-table(ifelse(subset_df[,2]>cut_off_methylation,"hyper","hypo"))
  stats_meth<-rbind(stats_meth,status_meth)
  
}
names(stats_meth)<-c("hyper","hypo")
input_barplot<-data.frame(group=c("MUT","MUT","MUT","WT","WT","WT"),melt(cbind(list_samples,stats_meth)))

pdf("GlobalMethLevels.1000kb.Barcharts.pdf")

p<-ggplot(input_barplot, aes(fill=variable, y=value, x=list_samples)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_manual(values = c("hyper" = "red3", "hypo" = "green3"))

print(p)

input_barplot2<-input_barplot[,-2]
input_barplot2_agg<-aggregate(.~group+variable, input_barplot2, median)

p2<-ggplot(input_barplot2_agg, aes(fill=variable, y=value, x=group)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_manual(values = c("hyper" = "red3", "hypo" = "green3"))

print(p2)

dev.off()

#
# Barplot methylation at group level
#



# png("GlobalMethLevels.1000kb.Violin.png")
# p <- ggplot(df_meth_all, aes(x=samples, y=meth_signals,color=samples)) + 
#   geom_violin(trim=FALSE)+ ylab("M/(M+U)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# print(p)
# dev.off()

# png("GlobalMethLevels.1000kb.yesredo.nocontrol.Violin.png")
# df_meth_all2<-df_meth_all[df_meth_all[,1]%in%c("IAP89","IAP89B","IAPA","IwtAP2R4a","IwtAP6c_redo","IWTAPB"),]
# p <- ggplot(df_meth_all2, aes(x=samples, y=meth_signals,color=samples)) + 
#   geom_violin(trim=FALSE)+ ylab("M/(M+U)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# print(p)
# dev.off()


# png("GlobalMethLevels.ViolinPlot.1000kb.png")
# dp <- ggplot(df_meth_all, aes(x=samples, y=meth_signals, fill=samples)) + 
# 	geom_violin(trim=FALSE)+
# 	geom_boxplot(width=0.1, fill="white")+ ylab("M/(M+U)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# print(dp)
# dev.off()

