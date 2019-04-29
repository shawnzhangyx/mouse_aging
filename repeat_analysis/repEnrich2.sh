#python2.7 RepEnrich2/RepEnrich2_subset.py D0.bam 30 D0 --pairedend TRUE
# 
for day in 3 10 18; do 
python2.7 RepEnrich2/RepEnrich2_subset.py ../../data/snATAC/bowtie2_bam/frontal_cortex_${day}_months.X2000.sorted.bam 30 ../../analysis/repeats_RepEnrich2/bam_subset/frontal_cortext_${day}months --pairedend TRUE
done

# Run RepEnrich2
python2.7 RepEnrich2/RepEnrich2.py \
  ~/software/RepEnrich2/refs/mm10_repeatmasker.txt \
  ../../analysis/repeats_RepEnrich2/ \
  frontal_cortex_3m  \
  ~/software/RepEnrich2/refs/setup_folder_mm10 \
  ../../analysis/repeats_RepEnrich2/bam_subset/frontal_cortext_3months_multimap_R1.fastq \
  --fastqfile2 ../../analysis/repeats_RepEnrich2/bam_subset/frontal_cortext_3months_multimap_R2.fastq \
  ../../analysis/repeats_RepEnrich2/bam_subset/frontal_cortext_3months_unique.bam \
  --cpus 16 --pairedend TRUE


