# The scripts in thie folder should be run after the snapATAC anaysis. 

## combine all up and down regulated peaks during aging.
cat ../../analysis/snapATAC/*/age_diff_edgeR.snap/*.up.bed > ../../analysis/snapATAC/all_celltypes/allcelltypes.age_up.peaks.bed
cat ../../analysis/snapATAC/*/age_diff_edgeR.snap/*.down.bed > ../../analysis/snapATAC/all_celltypes/allcelltypes.age_down.peaks.bed

########## All celltype processing #######
mkdir ../../analysis/snapATAC/all_celltypes
cd ../../analysis/snapATAC/
bams=*/bam.cluster/*.sorted.bam
featureCounts -a ../../data/snATAC/peaks/all_tissue.merged.peaks.saf -o all_celltypes/all_celltypes.counts $bams -F SAF -T 20 -O
cd -
Rscript all_celltype.clustering.r

##merge all individual cell type peaks into merged peaks.
peaks=$(ls ../../analysis/snapATAC/*/peak.cluster/*replicated_peaks.narrowPeak)
cat $peaks |bedtools sort|bedtools merge > ../../analysis/snapATAC/all_celltypes/allcelltypes.cluster_peaks.merged.bed

# annotate all peaks.
cd ../../analysis/snapATAC/all_celltypes/
intersectBed -a allcelltypes.cluster_peaks.merged.bed \
  -b ../../../data/encode/peaks/H3K27ac_merged_1k.bed -c |\
  intersectBed -a - -b ../../../data/encode/peaks/H3K27me3_merged_1k.bed -c |\
  intersectBed -a - -b ../../../data/encode/peaks/H3K4me3_merged_1k.bed -c |\
  intersectBed -a - -b ../../../data/encode/peaks/H3K4me1_merged_1k.bed -c |\
  intersectBed -a - -b ../../../data/encode/peaks/H3K9me3_merged_1k.bed -c |\
  intersectBed -a - -b ../../../data/encode/peaks/CTCF/CTCF_merged_1k.bed -c \
  > allcelltypes.cluster_peaks.merged.annotated.txt
cd -
Rscript annotate_elements.all_celltypeCluster.r

Rscript plot.cell_num_reads.r

