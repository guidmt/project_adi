library(GenomicRanges)

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

df_ann<-unique(df_ann)[,c(1,2,3,11)]


###################################################################################################################################################################################################
# Gene Levels analysis
###################################################################################################################################################################################################

done = FALSE

if(done!=TRUE){
  
  library(data.table)
  
  setwd("/mnt/data/lab/gmt_data/data_brain/adi/CytosineReports")
  
  list_meth_file<-dir()
  
  list_samples<-vector(mode='list',length(list_meth_file))
  
  for(i in 1:length(list_meth_file)){
    
    print(i)
    
    print("reading data")
    meth_table_current_mouse<-data.frame(fread(file=list_meth_file[i]))
    colnames(meth_table_current_mouse)<-c("chromosome","position","strand","count_methylated","count_unmethylated","c_context","trinucleotide_context")
    
    beta_current_values<-data.frame()
    
    #select only the sites in which the gens were identified DEGs
    for(row_ann in 1:nrow(df_ann)){
      
      print(row_ann)
      
      chrom<-paste("chr",df_ann[row_ann,1],sep="")
      start<-df_ann[row_ann,2]
      end<-df_ann[row_ann,3]
      genes<-df_ann[row_ann,ncol(df_ann)]
      
      test<-meth_table_current_mouse[which(meth_table_current_mouse[,1]==chrom &
                                             meth_table_current_mouse[,2]>=start & meth_table_current_mouse[,2]<=end),]
      
      if(nrow(test)>=1){
        
        coordinates<-paste(range(test[,1]),collapse="-")
        
        beta_values<-sum(test$count_methylated)/(sum(test$count_methylated)+sum(test$count_unmethylated))
        
        df_beta<-data.frame(gene_symbol=genes,beta_values)
        
        beta_current_values<-rbind(beta_current_values,df_beta)
        
      }
      
    }
    
    
    list_samples[[i]]<-beta_current_values
    
  }
  
  
  res_meth<-do.call(cbind,list_samples)
  
  res_meth2<-res_meth[,grep(colnames(res_meth),pattern="beta_values")]
  colnames(res_meth2)<-sapply(strsplit(list_meth_file,split="\\."),"[[",1)
  
  res_meth3<-as.data.frame(unique(cbind(genes=as.character(res_meth[,1]),res_meth2)))
  
  remove_inf<-function(x){
    idx<-which(!is.finite(x))
    x[idx]<-0
    return(x)
  }
  
  res_meth3[,-1]<-sapply(res_meth3[,-1],remove_inf)
  
  colnames(res_meth3)[1]<-"gene_symbol"
  
  
  setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")
  
  save(res_meth3,file="DynamicCall_methylation_data_DMET_CpG_30_03_2020_pval0.05_lfc1_newestimation.December.RData")
  
  count_check<-function(x){
    
    hypo<-length(x[x<=0.1])/length(x)
    partial_meth<-length(x[x>0.1 & x< 0.5])/length(x)
    hyper<-length(x[x>=0.5])/length(x)
    res<-data.frame(hypo,partial_meth,hyper)
  }
  
  pdf("Dynamic_beta_hist.pdf")
  
  for(i in 2:ncol(res_meth3[,-1])){
    
    test<-hist(res_meth3[,i],main=colnames(res_meth3)[i],breaks=30)
  }
  dev.off()
  
  pdf("DynamicCall_beta_hist_scatterplot.pdf")
  
  save_x<-vector(mode="list",length(2:ncol(res_meth3[,-1])))
  save_y<-vector(mode="list",length(2:ncol(res_meth3[,-1])))
  
  for(i in 2:ncol(res_meth3[,-1])){
    
    test<-hist(res_meth3[,i],main=colnames(res_meth3)[i],breaks=30)
    
    save_x[[i]]<-test$mids
    save_y[[i]]<-test$counts
  }
  names(save_x)<-colnames(res_meth3[,-1])
  names(save_y)<-colnames(res_meth3[,-1])
  
  save_x2<-save_x[c(2:5,7,8)]
  save_y2<-save_y[c(2:5,7,8)]
  
  
  colors<-c("red3","red2","red1","green3","green2","green1")
  
  pdf("DynamicCall_beta_hist_scatterplot.pdf")
  
  plot(save_x2[[1]],save_y2[[1]],ylim=c(0,2000),col="red3",type="l") 
  
  for(p in 2:length(save_x2)){
    par(new=TRUE)
    plot(save_x2[[p]],save_y2[[p]],ylim=c(0,2000),xaxt="n",yaxt="n",xlab="",ylab="",col=colors[p],type="l") 
  }
  dev.off()
  
  min_meth<-apply(res_meth3[-1],2,min)[-c(1,6)]
  mu_meth<-apply(res_meth3[-1],2,mean)[-c(1,6)]
  max_meth<-apply(res_meth3[-1],2,max)[-c(1,6)]
  
  stats_methylation<-do.call(rbind,apply(res_meth3[-1],2,count_check))[-c(1,6),]
  
  
  pdf("Dynamic_beta_dmrs_hist.pdf")
  
  barplot(t(stats_methylation), main="mean % methylated genes (collapsed regions)",
          xlab="Sample", col=c("green3","orange3","red3"),
          legend = rownames(t(stats_methylation)))
  
  mtx<-t(stats_methylation)
  mtx2<-rbind(mtx[1,],apply(mtx[2:3,],2,sum))
  rownames(mtx2)<-c("hypo","hyper")
  
  barplot(mtx2, main="mean % methylated genes (collapsed regions)",
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
  
  barplot(t(mtx3), main="mean % methylated genes (collapsed regions)",
          xlab="Sample", col=c("green3","red3"),
          ylab="mean % methylated genes (collapsed regions)",
          legend = rownames(t(mtx3)))
  
  dev.off()
  
}

###################################################################################################################################################################################################
# Region Levels analysis
###################################################################################################################################################################################################
library(GenomicRanges)

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

df_ann<-unique(df_ann)[,c(1,2,3,11)]


library(data.table)

setwd("/mnt/data/lab/gmt_data/data_brain/adi/CytosineReports")

list_meth_file<-dir()

list_samples<-vector(mode='list',length(list_meth_file))

for(i in 1:length(list_meth_file)){
  
  print(i)
  
  print("reading data")
  meth_table_current_mouse<-data.frame(fread(file=list_meth_file[i]))
  colnames(meth_table_current_mouse)<-c("chromosome","position","strand","count_methylated","count_unmethylated","c_context","trinucleotide_context")
  
  beta_current_values<-data.frame()
  
  #select only the sites in which the gens were identified DEGs
  for(row_ann in 1:nrow(df_ann)){
    
    print(row_ann)
    
    chrom<-paste("chr",df_ann[row_ann,1],sep="")
    start<-df_ann[row_ann,2]
    end<-df_ann[row_ann,3]
    genes<-df_ann[row_ann,ncol(df_ann)]
    
    test<-meth_table_current_mouse[which(meth_table_current_mouse[,1]==chrom &
                                           meth_table_current_mouse[,2]>=start & meth_table_current_mouse[,2]<=end),]
    
    if(nrow(test)>=1){
      
      coordinates<-paste(range(test[,1]),collapse="-")
      
      beta_values<-(test$count_methylated)/((test$count_methylated)+(test$count_unmethylated))
      
      df_beta<-data.frame(gene_symbol=genes,beta_values)
      
      beta_current_values<-rbind(beta_current_values,df_beta)
      
    }
    
  }
  
  
  list_samples[[i]]<-beta_current_values
  
}


res_meth<-do.call(cbind,list_samples)

res_meth2<-res_meth[,grep(colnames(res_meth),pattern="beta_values")]
colnames(res_meth2)<-sapply(strsplit(list_meth_file,split="\\."),"[[",1)

res_meth3<-as.data.frame(unique(cbind(genes=as.character(res_meth[,1]),res_meth2)))

remove_inf<-function(x){
  idx<-which(!is.finite(x))
  x[idx]<-0
  return(x)
}

res_meth3[,-1]<-sapply(res_meth3[,-1],remove_inf)

colnames(res_meth3)[1]<-"gene_symbol"


setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")

save(res_meth3,file="DynamicCall_methylation_data_DMET_CpG_30_03_2020_pval0.05_lfc1_newestimation.December.Regions.RData")

count_check<-function(x){
  
  hypo<-length(x[x<=0.1])/length(x)
  partial_meth<-length(x[x>0.1 & x< 0.5])/length(x)
  hyper<-length(x[x>=0.5])/length(x)
  res<-data.frame(hypo,partial_meth,hyper)
}

pdf("DynamicCall_beta_hist.regions.pdf")

for(i in 2:ncol(res_meth3[,-1])){
  
  test<-hist(res_meth3[,i],main=colnames(res_meth3)[i],breaks=30)
}
dev.off()


stats_methylation<-do.call(rbind,apply(res_meth3[-1],2,count_check))[-c(1,6),]


pdf("DynamicCall_beta_dmrs_hist.regions.pdf")

barplot(t(stats_methylation), main="mean % methylated regions",
        xlab="Sample", col=c("green3","orange3","red3"),
        legend = rownames(t(stats_methylation)))

mtx<-t(stats_methylation)
mtx2<-rbind(mtx[1,],apply(mtx[2:3,],2,sum))
rownames(mtx2)<-c("hypo","hyper")

barplot(mtx2, main="mean % methylated regions",
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
        ylab="mean % methylated regions",
        legend = rownames(t(mtx3)))

dev.off()
