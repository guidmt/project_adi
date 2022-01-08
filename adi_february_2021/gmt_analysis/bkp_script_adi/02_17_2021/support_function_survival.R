plotSurvival<-function(TCGA_surv,status,title){

	daysToDeath <- as.numeric(as.character(TCGA_surv[, "days_to_death"]))
	daysToLastFollowup <- as.numeric(as.character(TCGA_surv[, "days_to_last_follow_up"]))
	vitalStatus <- as.character(TCGA_surv[, "vital_status"])

	print(summary(daysToDeath - daysToLastFollowup))

	daysToEvent <- rep(NA, nrow(TCGA_surv))
	daysToEvent[which(vitalStatus %in% "Alive")] <- daysToLastFollowup[which(vitalStatus %in% "Alive")]
	daysToEvent[which(vitalStatus %in% "Dead")] <- daysToDeath[which(vitalStatus %in% "Dead")]

	eventStatus <- rep(NA, nrow(TCGA_surv))
	eventStatus[which(vitalStatus %in% "Alive")] <- 0
	eventStatus[which(vitalStatus %in% "Dead")] <- 1

	surv_final<-data.frame(daysToEvent,eventStatus,status=TCGA_surv[,which(colnames(TCGA_surv)%in%status)])

	print(head(surv_final))

	library(survival)
	library(survminer)

	sfit <- survminer::surv_fit(Surv(surv_final[,"daysToEvent"],surv_final[,"eventStatus"]) ~ status,data=surv_final)
	
	print("done")

	p<-ggsurvplot(sfit, data = surv_final,conf.int=FALSE, pval=TRUE, risk.table=TRUE,risk.table.height=0.3,title=title)

	print(p)
	
}



plotSurvival2<-function(TCGA_surv,score_to_use,suffix_output){
  
  daysToDeath <- as.numeric(as.character(TCGA_surv[, "days_to_death"]))
  daysToLastFollowup <- as.numeric(as.character(TCGA_surv[, "days_to_last_follow_up"]))
  vitalStatus <- as.character(TCGA_surv[, "vital_status"])
  
  print(summary(daysToDeath - daysToLastFollowup))
  
  daysToEvent <- rep(NA, nrow(TCGA_surv))
  daysToEvent[which(vitalStatus %in% "Alive")] <- daysToLastFollowup[which(vitalStatus %in% "Alive")]
  daysToEvent[which(vitalStatus %in% "Dead")] <- daysToDeath[which(vitalStatus %in% "Dead")]
  
  eventStatus <- rep(NA, nrow(TCGA_surv))
  eventStatus[which(vitalStatus %in% "Alive")] <- 0
  eventStatus[which(vitalStatus %in% "Dead")] <- 1
  
  score_to_use2<-score_to_use
    
  surv_final<-data.frame(daysToEvent,
                         eventStatus,
                         global_score=TCGA_surv[,which(colnames(TCGA_surv)%in%score_to_use2)])
  
  colnames(surv_final)[ncol(surv_final)]<-"score"
  
  res.cut <- surv_cutpoint(surv_final, time = "daysToEvent", event = "eventStatus",
                           variables = c("score"))
  
  p1<-plot(res.cut, palette = "npg")
  print(p1)
  value_cutoff<-round(res.cut$cutpoint[2],2)

  res.cat <- surv_categorize(res.cut)
    
  fit <- survfit(Surv(daysToEvent, eventStatus)~score, data = res.cat)

  title<-paste(suffix_output,score_to_use2,"CP =",round(res.cut$cutpoint[1],2),"S =",value_cutoff)
  
  p2<-ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE,pval=TRUE,newpage = FALSE,title=title)
  
  print(p2)
  
  res_new_category<-data.frame(TCGA_surv[,1],res.cat)
  
  return(res_new_category)[1]<-"sample_id"
}


coxphAnalysis<-function(TCGA_surv_coxph,mtx_exp,file_gsva,thr_cox=0.05){
  
daysToDeath <- as.numeric(as.character(TCGA_surv_coxph[, "days_to_death"]))
daysToLastFollowup <- as.numeric(as.character(TCGA_surv_coxph[, "days_to_last_follow_up"]))
vitalStatus <- as.character(TCGA_surv_coxph[, "vital_status"])
  
print(summary(daysToDeath - daysToLastFollowup))
  
daysToEvent <- rep(NA, nrow(TCGA_surv_coxph))
daysToEvent[which(vitalStatus %in% "Alive")] <- daysToLastFollowup[which(vitalStatus %in% "Alive")]
daysToEvent[which(vitalStatus %in% "Dead")] <- daysToDeath[which(vitalStatus %in% "Dead")]
  
eventStatus <- rep(NA, nrow(TCGA_surv_coxph))
eventStatus[which(vitalStatus %in% "Alive")] <- 0
eventStatus[which(vitalStatus %in% "Dead")] <- 1
  
surv_final_coxph<-data.frame(TCGA_surv_coxph[,1],daysToEvent,eventStatus)
surv_final_coxph[,1]<-as.character(surv_final_coxph[,1])
colnames(surv_final_coxph)[1]<-"ID"
  
mtx_genes_to_filter_coxph<-read.delim(file=file_gsva)

genes_to_filter_mtx_coxph<-as.character(mtx_genes_to_filter_coxph[,1])
  
mtx_ge_coxph<-data.frame(ID=colnames(mtx_exp),t(mtx_exp))
                           
mtx_ge_coxph2<-mtx_ge_coxph[,which(colnames(mtx_ge_coxph)%in%c("ID",genes_to_filter_mtx_coxph))]
                           
input_coxph<-merge(surv_final_coxph,mtx_ge_coxph2,by="ID")

require(RegParallel)

res <- RegParallel(
        data = input_coxph,
        formula = 'Surv(daysToEvent, eventStatus) ~ [*]',
        FUN = function(formula, data)
          coxph(formula = formula,
                data = data,
                ties = 'breslow',
                singular.ok = TRUE),
        FUNtype = 'coxph',
        variables = colnames(input_coxph)[4:ncol(input_coxph)],
        blocksize = 200,
        cores = 2,
        nestedParallel = FALSE,
        conflevel = 95)

res <- res[order(res$LogRank, decreasing = FALSE),]
finalCoxGenes <- subset(res, LogRank <= thr_cox)

finalCoxGenes2<-merge(finalCoxGenes,mtx_genes_to_filter_coxph,by.x="Variable",by.y="hs")
  
return(finalCoxGenes2)

}

                           