# link files
link_files.sh
# count reads. 
bams=$(ls shared.endo/bam/*.bam)
featureCounts -a ../../data/snATAC/peaks/all_tissue.merged.peaks.saf -o shared.endo/shared_endo.counts $bams -F SAF -T 20 -O 

