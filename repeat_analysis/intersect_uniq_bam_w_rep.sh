#sample=DH_03_rep2
for sample in DH_03_rep2 DH_10_rep2 DH_18_rep2; do 
intersectBed -a \
  /projects/ps-renlab/yanxiao/software/RepEnrich2/refs/setup_folder_mm10/repnames.bed \
  -b ../../data/snATAC/bam.filter.sort/$sample/$sample.filter.csort.bam \
  -wa -wb |cut -f 1-4,8|sort -k5,5 |gzip > ../../analysis/repeat_analysis/tmp.$sample.filter.repeat.txt.gz
done
