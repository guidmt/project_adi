library(openxlsx)

setwd("/Users/guidantonio/Desktop/may_adi2021/14_05_2021")

df2 <- unique(read.xlsx(xlsxFile = "C1AC1B_integrate_TME_1_0.05_v2.xlsx", sheet = 1))
df2_unique<-unique(df2[,c(1:12,14,15)])
write.xlsx(df2_unique,"C1AC1B_integrate_TME_1_0.05_v2_unique.xlsx",sep="\t",row.names=F,quote=F)

TCGA-DU-5872 TCGA-DU-7304 TCGA-FG-5965 
2            2            3 
