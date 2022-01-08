libLinear_mouse_hsApril<-function(input_for_ml3,input_mm2,list_ann=list_ann,ML="liblinear",prefix_output=NA,type_tumours="LGG",suffix_output,dwn_sampling=T){
  
  rownames(input_for_ml3)<-make.names(input_for_ml3[,1],unique=T)
  
  list_probabilities_mouse<-vector(mode="list",length(list_ann))
  list_probabilities_hs<-vector(mode="list",length(list_ann))
  all_performances<-data.frame()
  all_performances2<-data.frame()
  all_confmatrix<-data.frame()
  
  nsamples_test<-NULL
  
  for(status_ml in 1:length(list_ann)){
    
    status_id<-grep(colnames(input_for_ml3),pattern=paste("^",list_ann[status_ml],"$",sep=""))
    value_status_id<-grep(colnames(input_for_ml3),pattern=paste("^",list_ann[status_ml],"$",sep=""),value=T)
    
    input_for_ml3_temp<-input_for_ml3[which(input_for_ml3$type%in%type_tumours),]
    
    input_for_ml<-input_for_ml3_temp[,c(status_id,21:ncol(input_for_ml3))]
    
    if(length(grep(input_for_ml[,1],pattern="no_ann"))>0){
      input_for_ml<-input_for_ml[-which(input_for_ml[,1]%in%"no_ann"),]
    }else{
      input_for_ml<-input_for_ml
    }
    
    if(value_status_id=="IDH_TOT" & dwn_sampling == T){
      
      min_prop<-min(table(input_for_ml[,1]))
      
      #understand which is the class to use conserved
      which_min<-which.min(table(input_for_ml[,1]))
      which_max<-which.max(table(input_for_ml[,1]))
      
      name_class_to_preserve<-names(table(input_for_ml[,1]))[which_min]
      name_class_to_downsample<-names(table(input_for_ml[,1]))[which_max]
      
      input_for_ml_conserved<-input_for_ml[input_for_ml[,1]%in%name_class_to_preserve,]
      input_for_ml_downsample<-input_for_ml[sample(which(input_for_ml[,1]%in%name_class_to_downsample),min_prop),]
      
      input_for_ml<-rbind(input_for_ml_conserved,input_for_ml_downsample)
      
    }
    
    trainIndex <- createDataPartition(input_for_ml[,1], p = .6, 
                                      list = FALSE, 
                                      times = 1)
    
    train_matrix<-as.matrix(input_for_ml[trainIndex,-1])
    train_matrix[is.na(train_matrix)]<-0
    
    y_train<-input_for_ml[trainIndex,1]
    
    test<-input_for_ml[-trainIndex,-1]
    test[is.na(test)]<-0
    y_test<-input_for_ml[-trainIndex,1]
    
    print(status_ml)
    print(dim(test))
    
    co=heuristicC(train_matrix)
    
    if(value_status_id=="IDH_TOT"|value_status_id=="IDH1_R132H"){
      
      y_train<- relevel(as.factor(y_train), ref = "0")
      
    }else{
      
      y_train<-y_train
      
    }
    
    if(ML=="liblinear"){

    model <- LiblineaR(data=train_matrix,target=y_train,cost=co,type=0)
    predictions<-LiblineaR:::predict.LiblineaR(model, test,proba=T)$predictions
    prob_to_use <- LiblineaR:::predict.LiblineaR(model, test,proba=T)$probabilities
    prob_to_use<-data.frame(ID=as.character(rownames(test)),prob_to_use)
    
    }

    if(ML=="RF"){
    
    ctrl <- trainControl(method="cv",
			 number=5,
			 classProbs=T,
			 savePredictions=T,
			 allowParallel=TRUE)
    
    print("train")
    
    # train_matrix<-data.frame(y_train,train_matrix)
    # colnames(train_matrix)[1]<-"status"
    
    model <- train(y=y_train,
		     x=train_matrix,
		     method="ranger",
		     trControl=ctrl)

   predictions <- predict(model, test)
   prob_to_use<-predict(model, test,type="prob")
   prob_to_use<-data.frame(ID=as.character(rownames(test)),prob_to_use)

   }

    conf_matrix<-table(predictions,as.factor(y_test))
    rownames(conf_matrix)<-paste(rownames(conf_matrix),"_pred",sep="")
    colnames(conf_matrix)<-paste(colnames(conf_matrix),"_real",sep="")
    confmatmelt<- melt(conf_matrix)
    colnames(confmatmelt)[2]<-"real"
    
    all_confmatrix<-rbind(all_confmatrix,data.frame(comparison=list_ann[status_ml],confmatmelt))
    
    if(value_status_id=="IDH_TOT"){
  
    overall_classification<-confusionMatrix(predictions,as.factor(y_test),positive="1")$overall
    CM<-confusionMatrix(predictions,as.factor(y_test),positive="1")
    
    }else{
      
    overall_classification<-confusionMatrix(predictions,as.factor(y_test))$overall
    CM<-confusionMatrix(predictions,as.factor(y_test))
    
    }
  
  
    #average the stats between groups 
    if(ncol(as.matrix(CM$byClass))>=2){
    statsCM<-apply(CM$byClass,2,mean)
    }else{
    statsCM<-CM$byClass
    }

    all_performances<-rbind(all_performances,overall_classification)
    all_performances2<-rbind(all_performances2,statsCM)

    classification_onto_mouse <- LiblineaR:::predict.LiblineaR(model,input_mm2[,-1],proba=T)$probabilities
    rownames(classification_onto_mouse)<-rownames(input_mm2)
    
    list_probabilities_mouse[[status_ml]]<-classification_onto_mouse
    
    dt_test<-data.frame(ID=rownames(test),test)
    colnames(dt_test)[1]<-"ID_prob"
    colnames(prob_to_use)[-1]<-paste("variable",colnames(prob_to_use)[-1],sep="_")
    
    dt_test2<-merge(prob_to_use,dt_test,by.x="ID",by.y="ID_prob")
    list_probabilities_hs[[status_ml]]<-dt_test2
    nsamples_test<-c(nsamples_test,nrow(test))
  }
  
  colnames(all_performances)<-c("Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull","AccuracyPValue","McnemarPValue")
  colnames(all_performances2)<-c("Sensitivity","Specificity","Pos_Pred_Value","Neg_Pred_Value","Precision","Recall","F1","Prevalence","Detection_Rate","Detection_Prevalence","Balanced_Accuracy")
  rownames(all_performances)<-list_ann
  rownames(all_performances2)<-list_ann

  matrix_prediction_mouse_melt<-melt(t(do.call(cbind,list_probabilities_mouse)))
  
  library(ggplot2)
  library(grid)
  library(rowr)
  library(gridExtra)
  library(plyr)
  
  matrix_prediction_hs_melt <-rbind.fill(list_probabilities_hs)
  string_order<-c("ID",grep(colnames(matrix_prediction_hs_melt),pattern="variable",value=T))
  
  mat1<-matrix_prediction_hs_melt[,colnames(matrix_prediction_hs_melt)%in%string_order]
  mat2<-matrix_prediction_hs_melt[,-which(colnames(matrix_prediction_hs_melt)%in%string_order)]
  
  mat_final<-cbind(mat1,mat2)

  res_liblinear<-vector(mode="list",6)
  res_liblinear[[1]]<-do.call(cbind,list_probabilities_mouse)
  res_liblinear[[2]]<-mat_final
  res_liblinear[[3]]<-nsamples_test
  res_liblinear[[4]]<-all_performances
  res_liblinear[[5]]<-all_confmatrix
  res_liblinear[[6]]<-all_performances2

  return(res_liblinear)
}
