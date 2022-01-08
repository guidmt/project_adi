# project_adi

gene_expression_scripts: The most important script is run_compare_signatures.GENERAL.R, used to created heatmaps, retrieve annotations etc.

methylation_scripts: scripts for the methylation analysis
- run_quality_plot_March2021: annotate DMRs, create volcano plot, prepare input for deeptools
- run_annotation_list_dem2_March2021.R: annotate DMRs, create volcano plot, prepare input for deeptools
- Methylation_barplot_mouse_human.R: prepare the TCGA methylation data with the DMRs genes in mouse
- Methylation_barplot_mouse_human.random.R: prepare the TCGA methylation adta with random genes in mouse
- TCGA_to_TURCAN_new.R: predict the Turcan dataset using the TCGA data and the DMRs genes in mouse
- tf_avg_profile.PromotersExons.R: plot avg methylation signals of TFs overlapping promoters and enhancers / and all DMRs
- liblinear_MM_to_HS_April_meth.R: use TCGA dataset to measure the potential of prediction of the DMRs in IDH+ and Methylation Clusters
- liblinear_MM_to_HS_April_meth.nodown.R: the same of before without downsampling
- liblinear_MM_to_HS_April_meth.Random.R: use TCGA dataset to measure the potential of prediction of the DMRs in IDH+ and Methylation Clusters with random genes
- liblinear_MM_to_HS_April_meth.RandomWithBoot.R: use TCGA dataset to measure the potential of prediction of the DMRs in IDH+ and Methylation Clusters with random genes with bootstrap
- liblinear_MM_to_HS_April_meth.RandomWithBoot.TrimmedHard.R:use TCGA dataset to measure the potential of prediction of the DMRs in IDH+ and Methylation Clusters with random genes with bootstrap and binarizing the methylation data
- liblinear_MM_to_HS_April_meth.RandomWithBoot.TrimmedHypo.R: :use TCGA dataset to measure the potential of prediction of the DMRs in IDH+ and Methylation Clusters with random genes with bootstrap and trimming the methylation data hypomethylated
- meth_turcan_analysis.R: create methylation Turcan PCA

processing_meth_April/
- processing_full_genome.R: sliding window approache 1000 kb and measure methylation levels
- processing_full_genome.GenomicFeatures.R: like processing_full_genome.R for promoters and enhancers
- processing_meth_Adi_April_Dynamic.R: processing methylation signal using the selected DMRs Dynamically selected regions
- processing_meth_Adi_April.R: processing methylation signal using the selected DMRs regions
- run_annotation_list_dem2_April2021.R: QCs on methylation data

old_analysis: this folder contains the results that have been generated  from the analysis
