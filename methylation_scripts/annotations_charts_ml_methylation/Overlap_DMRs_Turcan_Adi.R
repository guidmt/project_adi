#
# See the overlap between the DMRs regions in Turcan and Adi
#
library(VennDetail)

#
#Supplementary Table 1. Differentially methylated genes in mutant IDH-expressing human astrocytes
#

require("biomaRt")
#Warning-use the last version of Ensembl to convert mgi symbol to hg symbol - call useEnsembl to use a different version
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#
# Turcan used a cut-off of 1 of "Normalized fold change (mut/wt)" to discriminate the UP and down DMRs
#

tab_dmrs_turcan<-read.delim(file="Turcan_DMRs_SuppTable1.txt",stringsAsFactors=F,fill=T)
up_dmrs_turcan<-unique(tab_dmrs_turcan[tab_dmrs_turcan[,3] >= 1,2])
up_dmrs_turcan<-unique(up_dmrs_turcan[up_dmrs_turcan!=""])

# it is possible that there are up and down methylated regions at the same genes, the have different cg code

down_dmrs_turcan<-unique(tab_dmrs_turcan[tab_dmrs_turcan[,3] <= -1,2])
down_dmrs_turcan<-unique(down_dmrs_turcan[down_dmrs_turcan!=""])

intersect(up_dmrs_turcan,down_dmrs_turcan)

tab_dmrs_Adi<-read.delim(file="annotation_dmrs_adi_pval0.05_lfc1.txt")
up_drms_adi<-unique(tab_dmrs_Adi[tab_dmrs_Adi$status%in%"DMRs-up","annot.symbol"])
hs_dmrs_up_from_mm = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = up_drms_adi , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

down_dmrs_adi<-unique(tab_dmrs_Adi[tab_dmrs_Adi$status%in%"DMRs-down","annot.symbol"])
hs_dmrs_down_from_mm = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = down_dmrs_adi , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

# also in the case of adi there is a problem of the same genes in common between up and down (16)
intersect(up_drms_adi,down_dmrs_adi)

ven <- venndetail(list(hyper_turcan = up_dmrs_turcan, 
                       hypo_turcan = down_dmrs_turcan,
                       hyper_adi = hs_dmrs_up_from_mm[,2],
                       hypo_adi = hs_dmrs_down_from_mm[,2]))

pdf("Venn_DMRs_Adi_Turcan.pdf")
plot(ven)
plot(ven, type = "upset")
dev.off()

Venn_Adi_Turcan<-result(ven, wide = TRUE)

write.table(Venn_Adi_Turcan,file="Venn_DMRs_Adi_Turcan.txt",sep="\t",row.names=F,quote=F)

