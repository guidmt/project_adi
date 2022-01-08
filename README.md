# project_adi

gene_expression_scripts: The most important script is run_compare_signatures.GENERAL.R, used to created heatmaps, retrieve annotations etc.

- run_compare_signatures.GENERAL.R: script to create the heatmap C1/C2 and select the genes. Create also the input for tumourmap. Use sampleinfo_C1vsC2.txt file. Works also with C1A/C1B
- run_survival_analysis_final.R: define the C1A/C1B groups and create and HTML reports. It works with several groups (IDH-no-codel, TRIPLE MUT,LGG, LGG + GBM)

/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/analysis_04_02_2021/res_04_02_2021/adi_february_2021/Comp8_C1A_vs_C1B/MultiOmicsSurv
- SurvWithMacroTME_CbioAdi.R: create overall survivals with a grouping defined by Adi and C1A/C1B classes
- SurvWithMacroTME.R: overall survivals C1A/C1B and TME microenvironment 
- SurvWithTMELasso.R: select the most important TME with lasso and do overall survival, integration C1/C1B.
- SurvWithTME.R: another script in which we integrate C1A/C1B and TME


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
