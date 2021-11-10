mkdir ../../analysis/snapATAC/de_peaks.subsample
# count all the bam files. 
cd ../../analysis/snapATAC/
bams=*/bam.cluster_age_rep/*.sorted.bam
featureCounts -a ../../data/snATAC/peaks/all_tissue.merged.peaks.saf -o de_peaks.subsample/all_celltypes.all_sample.counts $bams -F SAF -T 20 -O
cd -
#Rscript subsample_clusters_0.75M.r
Rscript subsample_clusters_0.5M.r

cd ../../analysis/snapATAC/de_peaks.subsample
bams=bam.subsample.0.5M/*.bam
featureCounts -a ../../../data/snATAC/peaks/all_tissue.merged.peaks.saf -o select_celltypes.0.5M.subsample.counts $bams -F SAF -T 20 -O
cd -

Rscript diff_analysis_edger.0.5M.r

cd ../../analysis/snapATAC/de_peaks.subsample
bams=bam.subsample.grids/*.bam
featureCounts -a ../../../data/snATAC/peaks/all_tissue.merged.peaks.saf -o select_celltypes.grid.subsample.counts $bams -F SAF -T 20 -O
cd -
Rscript diff_analysis_edger.grids.r
