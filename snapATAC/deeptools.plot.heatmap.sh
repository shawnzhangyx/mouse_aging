cd ../../analysis/snapATAC/FC/
computeMatrix reference-point \
  -S bigWig.cluster_age_rep/FC.metacell_1.*.bw \
  -R age_diff_edgeR.snap/1.up.bed \
  --referencePoint center \
  -a 2000 -b 2000 \
  -p 8 \
  -out FC.1.age_up_peaks.tab.gz 

plotHeatmap \
 -m FC.1.age_up_peaks.tab.gz\
 -out FC.1.age_up_peaks_signal.pdf \
 --heatmapHeight 15  \
 --plotTitle 'ATAC signal' 

computeMatrix reference-point \
  -S bigWig.cluster_age_rep/FC.metacell_1.*.bw \
  -R age_diff_edgeR.snap/1.down.bed \
  --referencePoint center \
  -a 2000 -b 2000 \
  -p 8 \
  -out FC.1.age_down_peaks.tab.gz  

plotHeatmap \
 -m FC.1.age_down_peaks.tab.gz\
 -out FC.1.age_down_peaks_signal.pdf \
 --heatmapHeight 15  \
 --plotTitle 'ATAC signal'  


