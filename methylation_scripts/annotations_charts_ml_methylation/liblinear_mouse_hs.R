libLinear_mouse_hs<-function(input_for_ml3,input_mm2,list_ann=list_ann,prefix_output=NA,type_tumours="LGG",suffix_output){
  
  rownames(input_for_ml3)<-paste("custom",1:nrow(input_for_ml3),sep="_")
  
  list_probabilities_mouse<-vector(mode="list",length(list_ann))
  list_probabilities_hs<-vector(mode="list",length(list_ann))
  all_performances<-data.frame()
  
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
    
    if(value_status_id=="IDH1_status"){
      
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
    print(dim(input_for_ml))
    print(dim(train_matrix))
    print(dim(test))

    co=heuristicC(train_matrix)
   
    if(value_status_id=="IDH1_status"|value_status_id=="IDH1_R132H"){

    y_train<- relevel(as.factor(y_train), ref = "0")
    
    }else{
    
    y_train<-y_train

    }


    model <- LiblineaR(data=train_matrix,target=y_train,cost=co,type=0)
    predictions<-LiblineaR:::predict.LiblineaR(model, test,proba=T)$predictions
    prob_to_use <- LiblineaR:::predict.LiblineaR(model, test,proba=T)$probabilities
    prob_to_use<-data.frame(ID=as.character(rownames(test)),prob_to_use)
      
    overall_classification<-confusionMatrix(predictions,as.factor(y_test))$overall
    all_performances<-rbind(all_performances,overall_classification)
    
    classification_onto_mouse <- LiblineaR:::predict.LiblineaR(model,input_mm2[,-1],proba=T)$probabilities
    rownames(classification_onto_mouse)<-rownames(input_mm2)
    
    list_probabilities_mouse[[status_ml]]<-classification_onto_mouse
    list_probabilities_hs[[status_ml]]<-prob_to_use
  }
  
  colnames(all_performances)<-c("Accuracy","Kappa","AccuracyLower","AccuracyUpper","AccuracyNull","AccuracyPValue","McnemarPValue")
  rownames(all_performances)<-list_ann
    
  matrix_prediction_mouse_melt<-melt(t(do.call(cbind,list_probabilities_mouse)))
  
  library(ggplot2)
  library(grid)
  library(gridExtra)
  
  pdf(paste(paste(prefix_output,"Hs_to_Mm",suffix_output,sep="_"),"pdf",sep="."))
  g1 <- ggplot(matrix_prediction_mouse_melt, aes(x = Var2, y=Var1)) +
    geom_tile(aes(fill = value))+scale_fill_gradient(
      low = "#FFFFFF",
      high = "#012345")
  print(g1)
  dev.off()
  
  library(rowr)

  matrix_prediction_hs_melt <- melt(Reduce(function(...) merge(..., by="ID", all=TRUE), list_probabilities_hs))
  matrix_prediction_hs_melt[is.na(matrix_prediction_hs_melt$value),3]<-0
  
  pdf(paste(paste(prefix_output,"Mm_to_Hs",suffix_output,sep="_"),"pdf",sep="."))
  g2 <- ggplot(matrix_prediction_hs_melt, aes(x = ID, y=variable)) +
    geom_tile(aes(fill = value))+scale_fill_gradient(
                                                     low = "#FFFFFF",
                                                     high = "#012345")
  print(g2)
  grid.newpage()
  grid.table(all_performances[,c(1,3,4)])
  dev.off()
  
  res_liblinear<-vector(mode="list",2)
  res_liblinear[[1]]<-do.call(cbind,list_probabilities_mouse)
  res_liblinear[[2]]<-do.call(cbind.fill,list_probabilities_hs)
  return(res_liblinear)
}
