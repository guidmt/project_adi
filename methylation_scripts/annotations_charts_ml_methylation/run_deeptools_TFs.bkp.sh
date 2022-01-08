#
# TP53 up and down 
#

computeMatrix reference-point --referencePoint center -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R TP53_regions_dmrs_adi_pval0.05_lfc1.REGIONS.txt -o TP53_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz \

plotHeatmap -m TP53_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out TP53_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf
plotProfile -m TP53_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out TP53_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf --yMin 0.21 --perGroup --legendLocation upper-right

#
# TP53 up  
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R TP53_regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.txt \
                                -o TP53_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz \

plotHeatmap -m TP53_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out TP53_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf
plotProfile -m TP53_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out TP53_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# TP53 down 
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R TP53_regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.txt \
                                -o TP53_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz \

plotHeatmap -m TP53_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out TP53_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf
plotProfile -m TP53_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out TP53_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right


# SUZ12 up and down 
#

computeMatrix reference-point --referencePoint center -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R SUZ12_regions_dmrs_adi_pval0.05_lfc1.REGIONS.txt -o SUZ12_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz \

plotHeatmap -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out SUZ12_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf
plotProfile -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out SUZ12_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf --yMin 0.23 --perGroup --legendLocation upper-right

#
# SUZ12 up  
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R SUZ12_regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.txt \
                                -o SUZ12_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz \

plotHeatmap -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out SUZ12_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf
plotProfile -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out SUZ12_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf --yMin 0.23 --perGroup --legendLocation upper-right

#
# SUZ12 down 
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R SUZ12_regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.txt \
                                -o SUZ12_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz \

plotHeatmap -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out SUZ12_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf
plotProfile -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out SUZ12_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# EZH2 up and down 
#

computeMatrix reference-point --referencePoint center -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R EZH2_regions_dmrs_adi_pval0.05_lfc1.REGIONS.txt -o EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz \

plotHeatmap -m EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf
plotProfile -m EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out EZH2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf  --yMin 0.22 --perGroup --legendLocation upper-right

#
# EZH2 up  
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R EZH2_regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.txt \
                                -o EZH2_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz \

plotHeatmap -m EZH2_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out EZH2_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf
plotProfile -m EZH2_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out EZH2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# EZH2 down 
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                  /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                  /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                  /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                  /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                  /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                  --beforeRegionStartLength 2000 \
                                  --afterRegionStartLength 2000 \
                                  --binSize 10 \
                                  -p 60 \
                                  --missingDataAsZero \
                                  -R EZH2_regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.txt \
                                  -o EZH2_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz \

plotHeatmap -m EZH2_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out EZH2_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf
plotProfile -m EZH2_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out EZH2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf --perGroup --legendLocation upper-right

#
# EOMES up and down 
#

computeMatrix reference-point --referencePoint center -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                              --beforeRegionStartLength 2000 \
                              --afterRegionStartLength 2000 \
                              --binSize 10 \
                              -p 60 \
                              --missingDataAsZero \
                              -R EOMES_regions_dmrs_adi_pval0.05_lfc1.REGIONS.txt -o EOMES_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz \

plotHeatmap -m EOMES_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out EOMES_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf
plotProfile -m EOMES_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out EOMES_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# EOMES up  
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R EOMES_regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.txt \
                                -o EOMES_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz \

plotHeatmap -m EOMES_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out EOMES_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf
plotProfile -m EOMES_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out EOMES_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# EOMES down 
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R EOMES_regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.txt \
                                -o EOMES_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz \

plotHeatmap -m EOMES_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out EOMES_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf
plotProfile -m EOMES_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out EOMES_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right


#
# FOXA1 up and down 
#

computeMatrix reference-point --referencePoint center -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                              --beforeRegionStartLength 2000 \
                              --afterRegionStartLength 2000 \
                              --binSize 10 \
                              -p 60 \
                              --missingDataAsZero \
                              -R FOXA1_regions_dmrs_adi_pval0.05_lfc1.REGIONS.txt -o FOXA1_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz \

plotHeatmap -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out FOXA1_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf
plotProfile -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out FOXA1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# FOXA1 up  
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R FOXA1_regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.txt \
                                -o FOXA1_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz \

plotHeatmap -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out FOXA1_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf
plotProfile -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out FOXA1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# FOXA1 down 
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R FOXA1_regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.txt \
                                -o FOXA1_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz \

plotHeatmap -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out FOXA1_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf
plotProfile -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out FOXA1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf --perGroup --legendLocation upper-right


#
# RUNX3 up and down 
#

computeMatrix reference-point --referencePoint center -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                              --beforeRegionStartLength 2000 \
                              --afterRegionStartLength 2000 \
                              --binSize 10 \
                              -p 60 \
                              --missingDataAsZero \
                              -R RUNX3_regions_dmrs_adi_pval0.05_lfc1.REGIONS.txt -o RUNX3_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz \

plotHeatmap -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out RUNX3_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf
plotProfile -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out RUNX3_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# RUNX3 up  
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R RUNX3_regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.txt \
                                -o RUNX3_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz \

plotHeatmap -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out RUNX3_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf
plotProfile -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out RUNX3_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# RUNX3 down 
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R RUNX3_regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.txt \
                                -o RUNX3_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz \

plotHeatmap -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out RUNX3_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf
plotProfile -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out RUNX3_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right


#
# SPIB up and down 
#

computeMatrix reference-point --referencePoint center -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                              --beforeRegionStartLength 2000 \
                              --afterRegionStartLength 2000 \
                              --binSize 10 \
                              -p 60 \
                              --missingDataAsZero \
                              -R SPIB_regions_dmrs_adi_pval0.05_lfc1.REGIONS.txt -o SPIB_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz \

plotHeatmap -m SPIB_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out SPIB_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf
plotProfile -m SPIB_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out SPIB_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# SPIB up  
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R SPIB_regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.txt \
                                -o SPIB_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz \

plotHeatmap -m SPIB_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out SPIB_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf
plotProfile -m SPIB_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out SPIB_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# SPIB down 
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R SPIB_regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.txt \
                                -o SPIB_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz \

plotHeatmap -m SPIB_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out SPIB_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf
plotProfile -m SPIB_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out SPIB_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf --perGroup --legendLocation upper-right


#
# TBX2 up and down 
#

computeMatrix reference-point --referencePoint center -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                              --beforeRegionStartLength 2000 \
                              --afterRegionStartLength 2000 \
                              --binSize 10 \
                              -p 60 \
                              --missingDataAsZero \
                              -R TBX2_regions_dmrs_adi_pval0.05_lfc1.REGIONS.txt -o TBX2_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz \

plotHeatmap -m TBX2_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out TBX2_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf
plotProfile -m TBX2_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out TBX2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# TBX2 up  
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R TBX2_regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.txt \
                                -o TBX2_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz \

plotHeatmap -m TBX2_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out TBX2_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf
plotProfile -m TBX2_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out TBX2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# TBX2 down 
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R TBX2_regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.txt \
                                -o TBX2_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz \

plotHeatmap -m TBX2_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out TBX2_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf
plotProfile -m TBX2_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out TBX2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right


#
# TFAP2A up and down 
#

computeMatrix reference-point --referencePoint center -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                              /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                              --beforeRegionStartLength 2000 \
                              --afterRegionStartLength 2000 \
                              --binSize 10 \
                              -p 60 \
                              --missingDataAsZero \
                              -R TFAP2A_regions_dmrs_adi_pval0.05_lfc1.REGIONS.txt -o TFAP2A_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz \

plotHeatmap -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out TFAP2A_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf
plotProfile -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out TFAP2A_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# TFAP2A up  
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R TFAP2A_regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.txt \
                                -o TFAP2A_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz \

plotHeatmap -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out TFAP2A_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf
plotProfile -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out TFAP2A_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# TFAP2A down 
#

computeMatrix reference-point --referencePoint center  -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
                                /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
                                --beforeRegionStartLength 2000 \
                                --afterRegionStartLength 2000 \
                                --binSize 10 \
                                -p 60 \
                                --missingDataAsZero \
                                -R TFAP2A_regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.txt \
                                -o TFAP2A_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz \

plotHeatmap -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out TFAP2A_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf
plotProfile -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out TFAP2A_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right



