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
                                -R TP53_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o TP53_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m TP53_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out TP53_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m TP53_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out TP53_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.21 --perGroup --legendLocation upper-right

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
                                -R TP53_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o TP53_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m TP53_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out TP53_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m TP53_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out TP53_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R TP53_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o TP53_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m TP53_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out TP53_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m TP53_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out TP53_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right


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
                                -R SUZ12_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o SUZ12_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out SUZ12_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out SUZ12_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.23 --perGroup --legendLocation upper-right

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
                                -R SUZ12_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o SUZ12_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out SUZ12_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out SUZ12_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.23 --perGroup --legendLocation upper-right

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
                                -R SUZ12_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o SUZ12_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out SUZ12_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m SUZ12_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out SUZ12_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

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
                                -R EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out EZH2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf  --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out EZH2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                  -R EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                  -o EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out EZH2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf --perGroup --legendLocation upper-right

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
                              -R EOMES_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o EOMES_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m EOMES_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out EOMES_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m EOMES_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out EOMES_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R EOMES_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o EOMES_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m EOMES_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out EOMES_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m EOMES_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out EOMES_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R EOMES_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o EOMES_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m EOMES_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out EOMES_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m EOMES_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out EOMES_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right


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
                              -R FOXA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o FOXA1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out FOXA1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out FOXA1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R FOXA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o FOXA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out FOXA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out FOXA1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R FOXA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o FOXA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out FOXA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m FOXA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out FOXA1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf --perGroup --legendLocation upper-right


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
                              -R RUNX3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o RUNX3_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out RUNX3_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out RUNX3_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R RUNX3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o RUNX3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out RUNX3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out RUNX3_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R RUNX3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o RUNX3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out RUNX3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m RUNX3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out RUNX3_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right


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
                              -R SPIB_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o SPIB_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m SPIB_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out SPIB_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m SPIB_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out SPIB_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R SPIB_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o SPIB_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m SPIB_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out SPIB_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m SPIB_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out SPIB_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R SPIB_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o SPIB_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m SPIB_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out SPIB_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m SPIB_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out SPIB_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf --perGroup --legendLocation upper-right


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
                              -R TBX2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o TBX2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m TBX2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out TBX2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m TBX2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out TBX2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R TBX2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o TBX2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m TBX2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out TBX2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m TBX2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out TBX2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R TBX2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o TBX2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m TBX2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out TBX2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m TBX2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out TBX2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right


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
                              -R TFAP2A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o TFAP2A_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out TFAP2A_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out TFAP2A_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R TFAP2A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o TFAP2A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out TFAP2A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out TFAP2A_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R TFAP2A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o TFAP2A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out TFAP2A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m TFAP2A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out TFAP2A_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# MYC up and down 
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
                              -R MYC_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o MYC_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m MYC_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out MYC_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m MYC_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out MYC_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# MYC up  
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
                                -R MYC_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o MYC_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m MYC_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out MYC_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m MYC_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out MYC_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# MYC down 
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
                                -R MYC_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o MYC_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m MYC_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out MYC_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m MYC_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out MYC_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# MYCN up and down 
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
                              -R MYCN_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o MYCN_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m MYCN_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out MYCN_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m MYCN_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out MYCN_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# MYCN up  
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
                                -R MYCN_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o MYCN_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m MYCN_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out MYCN_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m MYCN_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out MYCN_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# MYCN down 
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
                                -R MYCN_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o MYCN_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m MYCN_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out MYCN_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m MYCN_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out MYCN_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# E2F1 up and down 
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
                              -R E2F1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o E2F1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m E2F1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out E2F1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m E2F1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out E2F1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# E2F1 up  
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
                                -R E2F1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o E2F1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m E2F1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out E2F1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m E2F1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out E2F1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# E2F1 down 
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
                                -R E2F1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o E2F1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m E2F1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out E2F1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m E2F1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out E2F1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# E2F2 up and down 
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
                              -R E2F2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o E2F2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m E2F2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out E2F2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m E2F2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out E2F2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# E2F2 up  
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
                                -R E2F2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o E2F2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m E2F2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out E2F2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m E2F2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out E2F2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# E2F2 down 
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
                                -R E2F2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o E2F2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m E2F2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out E2F2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m E2F2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out E2F2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# E2F4 up and down 
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
                              -R E2F4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o E2F4_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m E2F4_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out E2F4_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m E2F4_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out E2F4_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# E2F4 up  
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
                                -R E2F4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o E2F4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m E2F4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out E2F4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m E2F4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out E2F4_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# E2F4 down 
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
                                -R E2F4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o E2F4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m E2F4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out E2F4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m E2F4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out E2F4_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# PRC2 up and down 
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
                              -R PRC2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o PRC2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m PRC2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out PRC2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m PRC2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out PRC2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# PRC2 up  
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
                                -R PRC2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o PRC2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m PRC2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out PRC2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m PRC2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out PRC2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# PRC2 down 
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
                                -R PRC2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o PRC2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m PRC2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out PRC2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m PRC2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out PRC2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# KDM4A up and down 
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
                              -R KDM4A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o KDM4A_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m KDM4A_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out KDM4A_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m KDM4A_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out KDM4A_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# KDM4AA up  
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
                                -R KDM4A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o KDM4A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m KDM4A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out KDM4A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m KDM4A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out KDM4A_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# KDM4A down 
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
                                -R KDM4A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o KDM4A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m KDM4A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out KDM4A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m KDM4A_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out KDM4A_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

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
                              -R EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m EZH2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out EZH2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# EZH2A up  
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
                                -R EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out EZH2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

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
                                -R EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m EZH2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out EZH2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right


#
# DMRT1 up and down 
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
                              -R DMRT1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o DMRT1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m DMRT1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out DMRT1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m DMRT1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out DMRT1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# DMRT1A up  
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
                                -R DMRT1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o DMRT1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m DMRT1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out DMRT1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m DMRT1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out DMRT1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# DMRT1 down 
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
                                -R DMRT1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o DMRT1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m DMRT1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out DMRT1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m DMRT1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out DMRT1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# GATA1 up and down 
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
                              -R GATA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o GATA1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m GATA1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out GATA1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m GATA1_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out GATA1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# GATA1 up  
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
                                -R GATA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o GATA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m GATA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out GATA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m GATA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out GATA1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# GATA1 down 
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
                                -R GATA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o GATA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m GATA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out GATA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m GATA1_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out GATA1_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# GATA2 up and down 
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
                              -R GATA2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o GATA2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m GATA2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out GATA2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m GATA2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out GATA2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# GATA2 up  
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
                                -R GATA2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o GATA2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m GATA2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out GATA2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m GATA2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out GATA2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# GATA2 down 
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
                                -R GATA2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o GATA2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m GATA2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out GATA2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m GATA2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out GATA2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# GATA3 up and down 
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
                              -R GATA3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o GATA3_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m GATA3_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out GATA3_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m GATA3_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out GATA3_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# GATA3 up  
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
                                -R GATA3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o GATA3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m GATA3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out GATA3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m GATA3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out GATA3_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# GATA3 down 
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
                                -R GATA3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o GATA3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m GATA3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out GATA3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m GATA3_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out GATA3_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right


#
# GATA4 up and down 
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
                              -R GATA4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o GATA4_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m GATA4_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out GATA4_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m GATA4_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out GATA4_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# GATA4 up  
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
                                -R GATA4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
                                -o GATA4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m GATA4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out GATA4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m GATA4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out GATA4_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# GATA4 down 
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
                                -R GATA4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
                                -o GATA4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m GATA4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out GATA4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m GATA4_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out GATA4_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# KDM2B up and down 
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
-R KDM2B_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o KDM2B_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m KDM2B_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out KDM2B_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m KDM2B_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out KDM2B_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# KDM2B up  
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
-R KDM2B_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
-o KDM2B_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m KDM2B_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out KDM2B_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m KDM2B_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out KDM2B_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# KDM2B down 
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
-R KDM2B_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
-o KDM2B_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m KDM2B_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out KDM2B_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m KDM2B_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out KDM2B_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

#
# OLIG2 up and down 
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
-R OLIG2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.txt -o OLIG2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz \

plotHeatmap -m OLIG2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out OLIG2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf
plotProfile -m OLIG2_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.gz -out OLIG2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.PromotersExons.center.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# OLIG2 up  
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
-R OLIG2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.UP.txt \
-o OLIG2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz \

plotHeatmap -m OLIG2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out OLIG2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf
plotProfile -m OLIG2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.gz -out OLIG2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.UP.withmiss.pdf --yMin 0.22 --perGroup --legendLocation upper-right

#
# OLIG2 down 
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
-R OLIG2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.DOWN.txt \
-o OLIG2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz \

plotHeatmap -m OLIG2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out OLIG2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf
plotProfile -m OLIG2_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.gz -out OLIG2_PlotProfile_regions_dmrs_adi_pval0.05_lfc1.PromotersExons.center.DOWN.withmiss.pdf  --perGroup --legendLocation upper-right

