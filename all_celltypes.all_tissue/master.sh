# call peaks from all tissue. 

# first count the reads 
for rank in {1..33}; do 
  ( echo $rank
  num=$(samtools view -c /projects/ps-renlab/lamaral/projects/Aging/all_tissues/split_bams/cluster_bams/$rank.metacell.bam)
  echo $rank $num >> ../../analysis/all_celltypes.all_tissue/cluster.read_counts.txt ) &
  done

Rscript select_qval_cut_offs.r

for i in {1..33}; do 
  line=$(sed -n ${i}p ../..//analysis/all_celltypes.all_tissue/cluster.read_counts.fdr.txt); 
  rank=$(echo $line|cut -f1 -d' ')
  fdr=$(echo $line|cut -f3 -d' ')
  echo $rank $fdr
  source activate py27 && macs2 callpeak -t /projects/ps-renlab/lamaral/projects/Aging/all_tissues/split_bams/cluster_bams/$rank.metacell.bam -n ../../analysis/all_celltypes.all_tissue/peaks.thresh/$rank -g mm --keep-dup all --shift 100 --extsize 200 --nomodel -q $fdr & 
  done 

# merge all peaks
cat ../../analysis/all_celltypes.all_tissue/peaks/*.narrowPeak |bedtools sort|mergeBed > ../../analysis/all_celltypes.all_tissue/all_celltypes.merged.peaks.bed
cat ../../analysis/all_celltypes.all_tissue/peaks.thresh/*.narrowPeak |bedtools sort|mergeBed > ../../analysis/all_celltypes.all_tissue/all_celltypes.merged.peaks.thresh.bed
intersectBed -a ../../analysis/all_celltypes.all_tissue/all_celltypes.merged.peaks.thresh.bed -b ../../annotations/mm10-blacklist.v2.bed -v |grep chrM -v > ../../analysis/all_celltypes.all_tissue/all_celltypes.merged.peaks.thresh.filter.bed


# overlap with histone marks and CTCF

# create a randomly generated set and overlap with the same features. 
cd ../../analysis/all_celltypes.all_tissue/
intersectBed -a all_celltypes.merged.peaks.thresh.bed \
  -b ../../data/encode/peaks/H3K27ac_merged_1k.bed -c |\
  intersectBed -a - -b ../../data/encode/peaks/H3K27me3_merged_1k.bed -c |\
  intersectBed -a - -b ../../data/encode/peaks/H3K4me3_merged_1k.bed -c |\
  intersectBed -a - -b ../../data/encode/peaks/H3K4me1_merged_1k.bed -c |\
  intersectBed -a - -b ../../data/encode/peaks/H3K9me3_merged_1k.bed -c |\
  intersectBed -a - -b ../../data/encode/peaks/CTCF/CTCF_merged_1k.bed -c \
  > allcelltypes.cluster_peaks.merged.thresh.annotated.txt

Rscript ~/software/github/seq-min-scripts/generate_random_segment_of_same_size.r \
  /projects/ps-renlab/share/bwa_indices/mm10.fa.fai \
  all_celltypes.merged.peaks.thresh.bed \
  control.peaks.thresh.bed 

intersectBed -a control.peaks.thresh.bed \
  -b ../../data/encode/peaks/H3K27ac_merged_1k.bed -c |\
  intersectBed -a - -b ../../data/encode/peaks/H3K27me3_merged_1k.bed -c |\
  intersectBed -a - -b ../../data/encode/peaks/H3K4me3_merged_1k.bed -c |\
  intersectBed -a - -b ../../data/encode/peaks/H3K4me1_merged_1k.bed -c |\
  intersectBed -a - -b ../../data/encode/peaks/H3K9me3_merged_1k.bed -c |\
  intersectBed -a - -b ../../data/encode/peaks/CTCF/CTCF_merged_1k.bed -c \
  > control.merged.thresh.annotated.txt

Rscript annotate_elements.all_celltypeCluster.r
cd -

bash ~/software/github/seq-min-scripts/bed_to_saf.sh ../../analysis/all_celltypes.all_tissue/all_celltypes.merged.peaks.thresh.filter.bed ../../analysis/all_celltypes.all_tissue/all_celltypes.merged.peaks.thresh.filter.saf
## count reads to these elements. 
cd /projects/ps-renlab/lamaral/projects/Aging/all_tissues/split_bams/cluster_bams/
bams=*.metacell.bam
featureCounts -a ~/projects/mouse_aging/analysis/all_celltypes.all_tissue/all_celltypes.merged.peaks.thresh.filter.saf -o ~/projects/mouse_aging/analysis/all_celltypes.all_tissue/all_celltypes.thresh.filter.counts $bams -F SAF -T 20 -O
cd -


