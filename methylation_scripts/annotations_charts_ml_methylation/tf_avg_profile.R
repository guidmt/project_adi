#
# Plot average methylation signals for the groups
#
library(profileplyr)
library(jsonlite)

tab_sampleinfo<-read.delim(file="/Users/guidantonio/Desktop/GitHub/project_adi/methylation_scripts/annotations_charts_ml_methylation/sampleinfo_TFs_plot.txt")

setwd("/Users/guidantonio/Desktop/GitHub/project_adi/methylation_scripts/annotations_charts_ml_methylation/avg_profile_inputs")

list_tf<-unique(tab_sampleinfo[,1])

# alphabetic order to have the DMRs as the first files to analyze
list_tf<-list_tf[order(list_tf)]

pdf("TF_avg_profiles.pdf",width=12,height=4)

for(plottf in list_tf){

 sub_ann_tf<- tab_sampleinfo[tab_sampleinfo[,1]%in%plottf,]

 nrow_ann<-nrow(sub_ann_tf)
 
 if(nrow_ann==3){
   
 sub_ann_tf<-sub_ann_tf[c(3,2,1),]
 
 }else{
   
 sub_ann_tf<-sub_ann_tf

 }
 

 par(mfrow=c(1,nrow_ann))
 
 for(exptf in 1:nrow(sub_ann_tf)){
   
  title_plot<- paste(plottf,sub_ann_tf[exptf,3])
  print(title_plot)
  
  proplyrObject <- import_deepToolsMat(as.character(sub_ann_tf[exptf,2]))
  
  myTempH <- scan(as.character(sub_ann_tf[exptf,2]), nlines = 1, what = "character", sep = "\n")
  myTemp <- read.delim(as.character(sub_ann_tf[exptf,2]), sep = "\t", comment.char = "@", 
                       header = FALSE)
  
  myTempM <- as.matrix(myTemp[, -seq(6)])
  if (sum(is.na(myTempM)) > 0) {
    myTempM[is.na(myTempM)] <- 0
    message("Matrix contained 'NA' values, these will be replaced with zeros")
  }
  
  myTempGR <- GRanges(myTemp[, 1], IRanges(myTemp[, 2], myTemp[, 
                                                               3]), names = myTemp[, 4], score = myTemp[, 5], strand = gsub("\\.", 
                                                                                                                            "*", myTemp[, 6]))
  
  
  info <- fromJSON(gsub("@", "", myTempH))
  bounds <- sample_boundaries <- info$sample_boundaries
  sample_labels <- info$sample_labels
  sample_boundariesToFilter <- info$sample_boundaries
  sample_Starts <- sample_boundaries[-length(sample_boundaries)]
  sample_Ends <- sample_boundaries[-1]
  sample_DataCol <- list()
  
  for (i in seq_along(sample_Starts)) {
    sample_DataCol[[i]] <- seq(sample_Starts[i] + 1, sample_Ends[i])
  }
  names(sample_DataCol) <- sample_labels
  
  #
  # This list contains the methylation profiles of each samples
  #
  myTempM_L <- lapply(sample_DataCol, function(x) myTempM[, 
                                                          x])
  myTempM_L <- lapply(myTempM_L, function(x) {
    colnames(x) <- NULL
    x
  })

  myTempM_L2<-lapply(myTempM_L,FUN=function(X){apply(X,2,mean)})
    
  myTempM_L3<-do.call(rbind,myTempM_L2)
  
  mut_samples_avg<-apply(myTempM_L3[1:3,],2,mean)
  wt_samples_avg<-apply(myTempM_L3[4:6,],2,mean)
  
  #plot(mut_samples_avg, type = 'l',col="red3",ylab="Methylation signal",xaxt="n",main=title_plot)
  #axis(1, at=c(0,100,200,300,400), labels=c("-2000","-1500","0","1500","2000"))
  #lines(wt_samples_avg, col = 'blue2')
  #legend("topright",  legend=c("MUT", "WT"), fill=c("red3","blue2"),col=c("red2","blue2"))
  

  mut_samples<-myTempM_L3[1:3,]
  wt_samples<-myTempM_L3[4:6,]

  matplot(t(mut_samples), type = 'l',col=c("orange3","salmon","orangered3"),ylab="Methylation signal",xaxt="n",main=title_plot,lty=c("dashed","dashed","dashed"))
  axis(1, at=c(0,100,200,300,400), labels=c("-2000","-1500","0","1500","2000"))
  lines(mut_samples_avg, col = 'darkred',lwd=2)
  matplot(t(wt_samples), type = 'l',col=c("royalblue4","lightblue4","steelblue1"),ylab="Methylation signal",xaxt="n",main=title_plot,add=T,lty=c("dashed","dashed","dashed"))
  lines(wt_samples_avg, col = 'darkblue',lwd=2)
  
  legend("topleft",   legend=c(colnames(t(mut_samples)),colnames(t(wt_samples))),fill=c("orange3","salmon","orangered3","royalblue4","lightblue4","steelblue1"),col=c("orange3","salmon","orangered3","royalblue4","lightblue4","steelblue1"))
  legend("topright",  legend=c("MUT", "WT"), fill=c("darkred","darkblue"),col=c("darkred","darkblue"))

  
 }
 
  
}
dev.off()

pdf("TF_avg_profiles.YScaleNew.pdf",width=12,height=4)

for(plottf in list_tf){
  
  sub_ann_tf<- tab_sampleinfo[tab_sampleinfo[,1]%in%plottf,]
  
  nrow_ann<-nrow(sub_ann_tf)
  
  if(nrow_ann==3){
    
    sub_ann_tf<-sub_ann_tf[c(3,2,1),]
    
  }else{
    
    sub_ann_tf<-sub_ann_tf
    
  }
  
  
  par(mfrow=c(1,nrow_ann))
  
  for(exptf in 1:nrow(sub_ann_tf)){
    
    title_plot<- paste(plottf,sub_ann_tf[exptf,3])
    print(title_plot)
    
    proplyrObject <- import_deepToolsMat(as.character(sub_ann_tf[exptf,2]))
    
    myTempH <- scan(as.character(sub_ann_tf[exptf,2]), nlines = 1, what = "character", sep = "\n")
    myTemp <- read.delim(as.character(sub_ann_tf[exptf,2]), sep = "\t", comment.char = "@", 
                         header = FALSE)
    
    myTempM <- as.matrix(myTemp[, -seq(6)])
    if (sum(is.na(myTempM)) > 0) {
      myTempM[is.na(myTempM)] <- 0
      message("Matrix contained 'NA' values, these will be replaced with zeros")
    }
    
    myTempGR <- GRanges(myTemp[, 1], IRanges(myTemp[, 2], myTemp[, 
                                                                 3]), names = myTemp[, 4], score = myTemp[, 5], strand = gsub("\\.", 
                                                                                                                              "*", myTemp[, 6]))
    
    
    info <- fromJSON(gsub("@", "", myTempH))
    bounds <- sample_boundaries <- info$sample_boundaries
    sample_labels <- info$sample_labels
    sample_boundariesToFilter <- info$sample_boundaries
    sample_Starts <- sample_boundaries[-length(sample_boundaries)]
    sample_Ends <- sample_boundaries[-1]
    sample_DataCol <- list()
    
    for (i in seq_along(sample_Starts)) {
      sample_DataCol[[i]] <- seq(sample_Starts[i] + 1, sample_Ends[i])
    }
    names(sample_DataCol) <- sample_labels
    
    #
    # This list contains the methylation profiles of each samples
    #
    myTempM_L <- lapply(sample_DataCol, function(x) myTempM[, 
                                                            x])
    myTempM_L <- lapply(myTempM_L, function(x) {
      colnames(x) <- NULL
      x
    })
    
    myTempM_L2<-lapply(myTempM_L,FUN=function(X){apply(X,2,mean)})
    
    myTempM_L3<-do.call(rbind,myTempM_L2)
    
    mut_samples_avg<-apply(myTempM_L3[1:3,],2,mean)
    wt_samples_avg<-apply(myTempM_L3[4:6,],2,mean)
    
    #plot(mut_samples_avg, type = 'l',col="red3",ylab="Methylation signal",xaxt="n",main=title_plot)
    #axis(1, at=c(0,100,200,300,400), labels=c("-2000","-1500","0","1500","2000"))
    #lines(wt_samples_avg, col = 'blue2')
    #legend("topright",  legend=c("MUT", "WT"), fill=c("red3","blue2"),col=c("red2","blue2"))
    
    
    mut_samples<-myTempM_L3[1:3,]
    wt_samples<-myTempM_L3[4:6,]
    
    matplot(t(mut_samples), type = 'l',col=c("orange3","salmon","orangered3"),ylab="Methylation signal",xaxt="n",main=title_plot,lty=c("dashed","dashed","dashed"),ylim=c(0.20,0.4))
    axis(1, at=c(0,100,200,300,400), labels=c("-2000","-1500","0","1500","2000"))
    lines(mut_samples_avg, col = 'darkred',lwd=2)
    matplot(t(wt_samples), type = 'l',col=c("royalblue4","lightblue4","steelblue1"),ylab="Methylation signal",xaxt="n",main=title_plot,add=T,lty=c("dashed","dashed","dashed"),ylim=c(0.20,0.4))
    lines(wt_samples_avg, col = 'darkblue',lwd=2)
    
    legend("topleft",   legend=c(colnames(t(mut_samples)),colnames(t(wt_samples))),fill=c("orange3","salmon","orangered3","royalblue4","lightblue4","steelblue1"),col=c("orange3","salmon","orangered3","royalblue4","lightblue4","steelblue1"))
    legend("topright",  legend=c("MUT", "WT"), fill=c("darkred","darkblue"),col=c("darkred","darkblue"))
    
    
  }
  
  
}
dev.off()
