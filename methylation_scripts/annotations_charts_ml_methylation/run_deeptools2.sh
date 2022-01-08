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
                              	-R regions_dmrs_adi_pval0.05_lfc1.REGIONS.bed -o regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz \

	plotHeatmap -m regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf
	plotProfile -m regions_dmrs_adi_pval0.05_lfc1.matrix.center.gz -out PlotProfile_regions_dmrs_adi_pval0.05_lfc1.matrix.center.withmiss.pdf --perGroup --legendLocation upper-right

#
# ALL UP with miss
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
                                -R regions_dmrs_adi_pval0.05_lfc1.REGIONS.UP.bed \
                                -o regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz \

plotHeatmap -m regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf
plotProfile -m regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.gz -out PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.UP.withmiss.pdf --perGroup --legendLocation upper-right

#
# ALL DOWN with miss
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
                              	-R regions_dmrs_adi_pval0.05_lfc1.REGIONS.DOWN.bed \
                              	-o regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz \
	
plotHeatmap -m regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf
plotProfile -m regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.gz -out PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.DOWN.withmiss.pdf --perGroup --legendLocation upper-right

#
# ALL DOWN PROMOTERS
#


#computeMatrix reference-point --referencePoint center -S /mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/ControlEC.bw \
#	/mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89B.bw \
#	/mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAP89.bw \
#	/mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IAPA.bw \
#	/mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP2R4a.bw \
#	/mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c.bw \
#	/mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IwtAP6c_redo.bw \
#	/mnt/data/lab/gmt_data/data_brain/Adi_Data-August28-2020/methVisualization/wigFiles/IWTAPB.bw \
#	--beforeRegionStartLength 1000 \
#	--afterRegionStartLength 1000 \
#	--binSize 10 \
#	-p 60 \
#	--missingDataAsZero \
#	-R regions_dmrs_adi_pval0.05_lfc1.REGIONS.Promoters.DOWN.bed \
#	-o regions_dmrs_adi_pval0.05_lfc1.center.Promoters.DOWN.gz \
#
#	plotHeatmap -m regions_dmrs_adi_pval0.05_lfc1.center.Promoters.DOWN.gz -out regions_dmrs_adi_pval0.05_lfc1.center.Promoters.DOWN.pdf
#plotProfile -m regions_dmrs_adi_pval0.05_lfc1.center.Promoters.DOWN.gz -out PlotProfile_regions_dmrs_adi_pval0.05_lfc1.center.Promoters.DOWN.pdf --perGroup --legendLocation upper-right


