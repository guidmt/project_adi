randomMatrixGenes<-function(matrix_hs = matrix_hs, matrix_all_genes){
  
  #the genes to use is equal to the number of DEGs genes
  number_genes_to_consider<-nrow(matrix_hs)
  #take the DEGs genes
  genes_degs<-rownames(matrix_hs)
  
  #remove the DEGs genes from the full matrix
  matrix_no_degs<-matrix_all_genes

  matrix_random_genes<-matrix_no_degs[sample(nrow(matrix_no_degs),number_genes_to_consider,replace=F),which(colnames(matrix_no_degs)%in%colnames(matrix_hs))]
  
  return(matrix_random_genes)
}


MLcore<-function(matrix_for_ml,annotation_for_classes,compute_conf_int=TRUE){
  
  #
  # ML1: Do  random forest
  #
  
  classes<-annotation_for_classes$status2
  classes<-ifelse(classes==1,"IDH_MUT","IDH_WT")
  
  input_ml<-data.frame(classes,t(matrix_for_ml))
  colnames(input_ml)[1]<-"status"
  input_ml$status<-factor(input_ml$status,levels=c("IDH_WT","IDH_MUT"))
  
  index <- createDataPartition(input_ml[,"status"], p = 0.7, list = FALSE)
  
  training<-input_ml[index,]
  test<-input_ml[-index,]
  
  print("trainControl")
  ctrl <- trainControl(method="cv",
                       number=5, 
                       classProbs=T, 
                       savePredictions=T,
                       sampling = "down",
                       summaryFunction=twoClassSummary,
                       allowParallel=TRUE)
  print("train")
  rfmodel <- train(status ~ .,
                   data=training,
                   method="ranger",
                   trControl=ctrl,
                   num.trees=100)
  

  rf_pred <- predict(rfmodel, test)
  
  print("measure performances")
  
  overall_per<-confusionMatrix(rf_pred, as.factor(test$status),positive="IDH_MUT")$overall
  byclass_per<-confusionMatrix(rf_pred, as.factor(test$status),positive="IDH_MUT")$byClass
  
  rf_pred <- predict(rfmodel, test,type="prob")
  
  rf.ROC <- roc(predictor=rf_pred$IDH_MUT,response=test$status)
  confidence_interval<-ci.auc(test$status, rf_pred$IDH_MUT, conf.level = 0.9) 
  ci_string<-paste(round(as.numeric(confidence_interval),3)[2],"(",round(as.numeric(confidence_interval),3)[1],"-",round(as.numeric(confidence_interval),3)[3],")",sep="")

  results_RF<-vector(mode="list",6)
  
  if(compute_conf_int==TRUE){
  
  ciobj <- ci.se(rf.ROC, specificities=seq(0, 1, l=25), conf.level=0.95, boot.n=10000)
  
  dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                       lower = ciobj[, 1],
                       upper = ciobj[, 3])
  
  results_RF[[1]]<-rfmodel
  results_RF[[2]]<-overall_per
  results_RF[[3]]<-byclass_per
  results_RF[[4]]<-rf.ROC
  results_RF[[5]]<-ci_string
  results_RF[[6]]<-dat.ci
  
  }else{
    
    results_RF[[1]]<-rfmodel
    results_RF[[2]]<-overall_per
    results_RF[[3]]<-byclass_per
    results_RF[[4]]<-rf.ROC
    results_RF[[5]]<-ci_string
    results_RF[[6]]<-"nothing"
    
  }
  
  return(results_RF)
  
}


permutationTestPvalue<-function(matrix_for_ml,annotation_for_classes,actual_auc,iter_perm){
  
  AUC_perm<-NULL
  
  for(i in 1:iter_perm){
  
  print(paste("permutation:",i))
    
  classes<-annotation_for_classes$status2
  classes<-ifelse(classes==1,"IDH_MUT","IDH_WT")
  
  input_ml<-data.frame(classes,t(matrix_for_ml))
  colnames(input_ml)[1]<-"status"
  
  input_ml$status<-sample(input_ml$status)
  
  input_ml$status<-factor(input_ml$status,levels=c("IDH_WT","IDH_MUT"))
  
  index <- createDataPartition(input_ml[,"status"], p = 0.7, list = FALSE)
  
  training<-input_ml[index,]
  test<-input_ml[-index,]
  
  print("trainControl")
  ctrl <- trainControl(method="cv",
                       number=5, 
                       classProbs=T, 
                       savePredictions=T,
                       sampling = "down",
                       summaryFunction=twoClassSummary,
                       allowParallel=TRUE)
  print("train")
  rfmodel <- train(biological_states ~ .,
                   data=training,
                   method="ranger",
                   trControl=ctrl,
                   num.trees=100)
  
  rf_pred <- predict(rfmodel, test)
  
  print("measure performances")
  
  overall_per<-confusionMatrix(rf_pred, as.factor(test$status),positive="IDH_MUT")$overall
  byclass_per<-confusionMatrix(rf_pred, as.factor(test$status),positive="IDH_MUT")$byClass
  
  rf_pred <- predict(rfmodel, test,type="prob")
  
  rf.ROC <- roc(predictor=rf_pred$IDH_MUT,response=test$status)
  
  auc_perm<-as.numeric(rf.ROC$auc)
  
  AUC_perm<-c(AUC_perm,auc_perm)
  
  }
  
  pvalue_auc_permutation<-length(which(AUC_perm>=actual_auc))/length(AUC_perm)
  
  res_permutation<-vector(mode="list",2)
  
  res_permutation[[1]]<-AUC_perm
  res_permutation[[2]]<-pvalue_auc_permutation
  
  return(res_permutation)
  
}


mlForModel<-function(matrix_hs,matrix_all_genes,annotation_for_classes,output_string,randGenes=F,permutation=F,iter_perm=100){
  
  #
  # ML1: Do  random forest
  #

  rfmodel<-MLcore(matrix_hs,annotation_for_classes)
  
  rf.ROC<-rfmodel[[4]]
  overall_per<-rfmodel[[2]]
  byclass_per<-rfmodel[[3]]
  ci_string<-rfmodel[[5]]
  dat.ci<-rfmodel[[6]]

  #
  # ML2: Do  a permutation to compute a p-value of the AUC
  #
  if(permutation==TRUE){

  actual_auc<-as.numeric(rf.ROC$auc)
  
  pvalueAUC<-permutationTestPvalue(matrix_hs,annotation_for_classes,iter_perm=iter_perm,actual_auc=actual_auc)
  
  output_hist<-paste(paste("histogram_auc_permutations",output_string,sep="."),"pdf",sep=".")
  
  pdf(output_hist)
  hist(pvalueAUC[[1]],xlim=c(0,1),main=round(pvalueAUC[[2]],3))
  abline(v=actual_auc,col="red",lty=3)
  dev.off()
  
  }else{
    
  pvalueAUC<-"Not computed"
    
  }
  
  output_roc<-paste(paste("ROC_randomforest",output_string,sep="."),"pdf",sep=".")
  
  pdf(output_roc)
  plot(rf.ROC,main=ci_string)
  require(gridExtra)
  per1<-tableGrob(data.frame(overall_per))
  per2<-tableGrob(data.frame(byclass_per))
  grid.arrange(per1,per2)
  
  pAUCCI<-ggroc(rf.ROC) + theme_minimal() + 
    geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + coord_equal() + 
    geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha= 0.2) + ggtitle(ci_string)
  
  print(pAUCCI)
  
  dev.off()
  
  #
  # ML1: Do  random forest random genes
  #
  
  if(randGenes==TRUE){
    
  output_roc<-paste(paste("ROC_randomforest",output_string,sep="."),".RANDOMGENES.pdf",sep=".")
  
  pdf(output_roc)
  
  print(plot(rf.ROC,col="grey"))

  AUC_average<-NULL
    
  for(i in 1:100){
    
  print(paste("Randomize genes stage:",i))
    
  matrix_random_genes_hs<-randomMatrixGenes(matrix_hs,matrix_all_genes)
    
  rfmodel_rand<-MLcore(matrix_random_genes_hs,annotation_for_classes,compute_conf_int=FALSE)
  
  rf.ROC_rand<-rfmodel_rand[[4]]
  overall_per_rand<-rfmodel_rand[[2]]
  byclass_per_rand<-rfmodel_rand[[3]]
  ci_string_rand<-rfmodel_rand[[5]]
  dat.ci_rand<-rfmodel_rand[[6]]
  
  print(plot(rf.ROC_rand,add=T,col="grey",type="l", lty=2))

  auc_current<-as.numeric(rf.ROC_rand$auc)
  
  AUC_average<-c(AUC_average,auc_current)
  
  }
  
  print(plot(rf.ROC,col="red",add=T,main=paste("AUC",round(mean(AUC_average)))))

  dev.off()
  
  }
  
  # output_roc<-paste(paste("ROC_randomforest",output_string,sep="."),"pdf",sep=".")
  # 
  # pdf(output_roc)
  # 
  # pAUCCI<-ggroc(rf.ROC_rand) + theme_minimal() + 
  #   geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + coord_equal() + 
  #   geom_ribbon(data = dat.ci_rand, aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha= 0.2) + ggtitle(ci_string_rand)
  # 
  # print(pAUCCI)
  # 
  # dev.off()
  
  
}
