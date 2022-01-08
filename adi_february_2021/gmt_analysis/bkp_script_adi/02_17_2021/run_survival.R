library(TCGAbiolinks)

setwd("/home/guidantoniomt/pseudospace/HMM")

tab_hmm<-read.delim(file="HMM_results_nstates_3.txt",stringsAsFactors=F)

list_tissue<-paste("TCGA-",unique(sapply(strsplit(tab_hmm[,1],split="\\."),"[[",1)),sep="")
tab_hmm[,1]<-unlist(lapply(strsplit(tab_hmm[,1],split="\\."),FUN=function(x){paste(x[2:4],collapse="-")}))


TCGA_global_tcga<-vector(mode="list",length(list_tissue))

for(tissue in 1:length(list_tissue)){

  print(tissue)

  clin.tissue<-GDCquery_clinic(list_tissue[tissue], "clinical")
  
  clin.tissue2<-clin.tissue[,which(colnames(clin.tissue)%in%c("submitter_id","vital_status", "days_to_death", "days_to_last_follow_up"))]

  # days_to_follow_up: Number of days between the date used for index and the date of the patient's last follow-up appointment or contact.
  # 

  TCGA_global_tcga[[tissue]]<-clin.tissue2

}


TCGA_global_tcga2<-do.call(rbind,TCGA_global_tcga)
TCGA_surv<-merge(TCGA_global_tcga2,tab_hmm[,c(1,5)],by.x="submitter_id",by.y="samples")

setwd("/home/guidantoniomt/pseudospace/survival_analysis")

pdf("test.pdf")
TCGAanalyze_survival(TCGA_surv,
                     "biological_states",
                     main = "TCGA",height = 10, width=10)
dev.off()

# source : https://bioinformatics.mdanderson.org/Supplements/ResidualDisease/Reports/assembleTCGAClinical.html

daysToDeath <- as.numeric(as.character(TCGA_surv[, "days_to_death"]))
daysToLastFollowup <- as.numeric(as.character(TCGA_surv[, "days_to_last_follow_up"]))
vitalStatus <- as.character(TCGA_surv[, "vital_status"])

summary(daysToDeath - daysToLastFollowup)

daysToEvent <- rep(NA, nrow(TCGA_surv))
daysToEvent[which(vitalStatus %in% "Alive")] <- daysToLastFollowup[which(vitalStatus %in% "Alive")]
daysToEvent[which(vitalStatus %in% "Dead")] <- daysToDeath[which(vitalStatus %in% "Dead")]

eventStatus <- rep(NA, nrow(TCGA_surv))
eventStatus[which(vitalStatus %in% "Alive")] <- 0
eventStatus[which(vitalStatus %in% "Dead")] <- 1

surv_final<-data.frame(daysToEvent,eventStatus,TCGA_surv$biological_states)

colnames(surv_final)[3]<-"biological_states"

library(survival)
library(survminer)

sfit <- survfit(Surv(as.numeric(surv_final$daysToEvent)/365, event = surv_final$eventStatus) ~ surv_final$biological_states,data=surv_final)

setwd("/home/guidantoniomt/pseudospace/survival_analysis")

pdf("HMM_nstates3_death.check_labels.pdf")
p1<-ggsurvplot(sfit)
print(p1)
dev.off()    

pdf("HMM_nstates3_death.pdf")
p2<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
           legend.labs=c("epi", "mes","mes-altered"), legend.title="biological_states",
           palette=c("dodgerblue2", "red2","orange2"),
           title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
           risk.table.height=0.3)
print(p2)
dev.off()

surv_final<-data.frame(daysToEvent,eventStatus,as.character(TCGA_surv$biological_states))
surv_final[,3]<-as.character(surv_final[,3])
surv_final[which(surv_final$biological_states%in%"mes"),3]<-"all_mes"
surv_final[which(surv_final$biological_states%in%"mix"),3]<-"all_mes"

sfit <- survfit(Surv(as.numeric(surv_final$daysToEvent)/365, event = surv_final$eventStatus) ~ surv_final$biological_states,data=surv_final)

pdf("HMM_nstates3_death_allmes_vs_epi.check_labels.pdf")
p1<-ggsurvplot(sfit)
print(p1)
dev.off()

pdf("HMM_nstates3_death_allmes_vs_epi.pdf")
p2<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
           legend.labs=c("all_mes_like","epi"), legend.title="biological_states",
           palette=c("red2","dodgerblue2"),
           title="Kaplan-Meier Curve for Epi, and Mes + Mes-altered",
           risk.table.height=0.3)
print(p2)
dev.off()


#
# Trim more than 5 years
#

surv_final<-data.frame(daysToEvent,eventStatus,TCGA_surv$biological_states)
colnames(surv_final)[3]<-"biological_states"
surv_final_trim<-surv_final[which(surv_final$daysToEvent/365 <=5),]

sfit <- survfit(Surv(as.numeric(surv_final_trim$daysToEvent)/365, event = surv_final_trim$eventStatus) ~ surv_final_trim$biological_states,data=surv_final_trim)

setwd("/home/guidantoniomt/pseudospace/survival_analysis")

pdf("HMM_nstates3_death.check_labels_trim5.pdf")
p3<-ggsurvplot(sfit)
print(p3)
dev.off()

pdf("HMM_nstates3_death_trim5.pdf")
p4<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
           legend.labs=c("epi", "mes","mes-altered"), legend.title="biological_states",
           palette=c("dodgerblue2", "red2","orange2"),
           title="Kaplan-Meier Curve for Epi, Mes, Mes-altered",
           risk.table.height=0.3)
print(p4)
dev.off()


surv_final<-data.frame(daysToEvent,eventStatus,as.character(TCGA_surv$biological_states))
surv_final[,3]<-as.character(surv_final[,3])
colnames(surv_final)[3]<-"biological_states"
surv_final_trim<-surv_final[which(surv_final$daysToEvent/365 <=5),]
surv_final_trim[which(surv_final_trim$biological_states%in%"mes"),3]<-"all_mes"
surv_final_trim[which(surv_final_trim$biological_states%in%"mix"),3]<-"all_mes"

sfit <- survfit(Surv(as.numeric(surv_final_trim$daysToEvent)/365, event = surv_final_trim$eventStatus) ~ surv_final_trim$biological_states,data=surv_final_trim)

pdf("HMM_nstates3_death_allmes_vs_epi.trim5.check_labels.pdf")
p1<-ggsurvplot(sfit)
print(p1)
dev.off()

pdf("HMM_nstates3_death_allmes_vs_epi.trim5.pdf")
p2<-ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE,
           legend.labs=c("all_mes_like","epi"), legend.title="biological_states",
           palette=c("red2","dodgerblue2"),
           title="Kaplan-Meier Curve for Epi, and Mes + Mes-altered",
           risk.table.height=0.3)
print(p2)
dev.off()

