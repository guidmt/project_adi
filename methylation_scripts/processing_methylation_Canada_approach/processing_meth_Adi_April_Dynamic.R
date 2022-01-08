library(GenomicRanges)
library(data.table)
library(tidyr)

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis")
source("support_function_heatmap_mouse.ge.R")

setwd("/mnt/data/lab/gmt_data/data_brain/adi/BriefSummaryReports")

input_met_mm<-read.delim(file="ListOfDifferentiallyMethylatedCGI-Mutants.txt",header=F)
input_met_mm_rid<-input_met_mm[,c(1:3)]
colnames(input_met_mm_rid)<-c("chr","start","end")
input_met_mm_rid[,1]<-gsub(input_met_mm_rid[,1],pattern="chr",replacement="")

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis")
annotation<-read.delim(file="RefSeq_mm_38_mm10.txt",header=T)
annotation2<-annotation[grep(annotation[,2],pattern="NM"),]
annotation2[,3]<-gsub(annotation2[,3],pattern="chr",replacement="")

annotation3<-annotation2[,c(13,3,5,6)]
colnames(annotation3)<-c("gene_symbol","chr","start","end")

meth_data<-makeGRangesFromDataFrame(input_met_mm_rid)
ann_data<-makeGRangesFromDataFrame(annotation3,keep.extra.columns = TRUE)

annotation<-findOverlaps(ann_data,meth_data,maxgap=1000)

qh<-queryHits(annotation)
sh<-subjectHits(annotation)

res_subject<-data.frame(meth_data[sh,])
res_query<-data.frame(ann_data[qh,])

df_ann <- cbind(res_subject,res_query)
colnames(df_ann)[2:3]<-c("Cpg_start","Cpg_end")

df_ann<-unique(df_ann)[,c(1,2,3,4)]

setwd("/mnt/data/lab/gmt_data/data_brain/adi/CytosineReports")

list_meth_file<-dir()

#Read each file
list_samples<-vector(mode='list',length(list_meth_file))

for(i in 1:length(list_meth_file)){
  
  print(i)
  
  print("reading data")
  meth_table_current_mouse<-data.frame(fread(file=list_meth_file[i]))
  colnames(meth_table_current_mouse)<-c("chromosome","position","strand","count_methylated","count_unmethylated","c_context","trinucleotide_context")
  
  # Step1: For the read.bismark function of bsseq, the cov is given by Cov[subjectHits(ol)] <- dt[queryHits(ol), list(Cov = (M + U))][["Cov"]]
  cov <- meth_table_current_mouse$count_methylated +  meth_table_current_mouse$count_unmethylated
  meth_table_current_mouse<-cbind(meth_table_current_mouse,cov)
  
  # Step2: the methylation level at a CpG site was the number of reads with that site methylated divided by the total number of reads covering the site.
  meth_signal <- (meth_table_current_mouse$count_methylated)/(meth_table_current_mouse$count_methylated+meth_table_current_mouse$count_unmethylated)
  
  meth_table_current_mouse<-cbind(meth_table_current_mouse,meth_signal)
  
  #only CpGs covered by ≥ 5x in all samples were retained for the computation of DNA methylation levels.
  idx_trim<- which(meth_table_current_mouse$cov > 5)
  
  meth_table_current_mouse2<-meth_table_current_mouse[idx_trim, ]
  
  
  #select only the sites in which the genes were identified DEGs
  wmeth_current_values<-data.frame()
  
  for(row_ann in 1:nrow(df_ann)){
    
    chrom<-paste("chr",df_ann[row_ann,1],sep="")
    start<-df_ann[row_ann,2]
    end<-df_ann[row_ann,3]
    genes<-df_ann[row_ann,ncol(df_ann)]
    
    # average methylation for a genomic region as the coverage-weighted 
    idx2<-which(meth_table_current_mouse2[,1]==chrom & meth_table_current_mouse2[,2]>=start & meth_table_current_mouse2[,2]<=end)
    
    test<-meth_table_current_mouse2[idx2,]
    
    # Step2: average methylation for a genomic region as the coverage-weighted 
    
    if(nrow(test)>=1){
      
      coordinates<-paste(range(test[,1]),collapse="-")
      
      wmeth<-weighted.mean(x=test$meth_signal,w=test$cov)
      
      df_wmeth<-data.frame(samples=list_meth_file[i],gene_symbol=genes,wmeth)
      wmeth_current_values<-rbind(wmeth_current_values,df_wmeth)
      
    }
    
  }
  
  
  list_samples[[i]]<-wmeth_current_values
  
}

res_meth<-do.call(rbind,list_samples)

res_meth2<-aggregate(.~samples+gene_symbol,res_meth ,mean)
res_meth_wide<-spread(res_meth2, key = samples, value = wmeth)
colnames(res_meth_wide)<-sapply(strsplit(colnames(res_meth_wide),split="\\."),"[[",1)

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/processing_meth_April")

save(res_meth_wide,file="methylation_levels_DMRs_Dynamic_approach_proc_Canada.RData")

pdf("weighthed_mean_DMRs_Dynamic_approach.pdf")
for(i in 2:ncol(res_meth_wide[,-1])){
  
  test<-hist(res_meth_wide[,i],main=colnames(res_meth_wide)[i],breaks=30)
}
dev.off()

count_check<-function(x){
  
  x<-x[!is.na(x)]
  hypo<-length(x[x<=0.1])/length(x)
  partial_meth<-length(x[x>0.1 & x< 0.5])/length(x)
  hyper<-length(x[x>=0.5])/length(x)
  res<-data.frame(hypo,partial_meth,hyper)
  
}

stats_methylation<-do.call(rbind,apply(res_meth_wide[-1],2,count_check))[-c(1,6),]

pdf("weighthed_mean_dmrs_hist.Dynamic_approach.regions.pdf")

barplot(t(stats_methylation), main="mean % methylated regions",
        xlab="Sample", col=c("green3","orange3","red3"),
        legend = rownames(t(stats_methylation)))

mtx<-t(stats_methylation)
mtx2<-rbind(mtx[1,],apply(mtx[2:3,],2,sum))
rownames(mtx2)<-c("hypo","hyper")

barplot(mtx2, main="mean % methylated genes regions",
        xlab="Sample", col=c("green3","red3"),
        legend = rownames(mtx2))

hypo_mut<-mean(mtx2[1,1:3])
hypo_wt<-mean(mtx2[1,4:6])

hyper_mut<-mean(mtx2[2,1:3])
hyper_wt<-mean(mtx2[2,4:6])

wt<-c(hypo_wt,hyper_wt)
mut<-c(hypo_mut,hyper_mut)

mtx3<-rbind(wt,mut)
colnames(mtx3)<-c("hypo","hyper")

barplot(t(mtx3), main="mean % methylated regions",
        xlab="Sample", col=c("green3","red3"),
        ylab="mean % regions",
        legend = rownames(t(mtx3)))

dev.off()

