tissue=$1
#sample=DH_03_rep2

for age in 03 10 18;do
  for rep in rep1 rep2; do
    sample=${tissue}_${age}_${rep}
    intersectBed -a \
  /projects/ps-renlab/yanxiao/software/RepEnrich2/refs/RepEnrich2_setup_mm10/repnames.bed \
  -b ../../data/snATAC/bam.filter.sort/$sample/$sample.filter.csort.bam \
  -wa -wb |cut -f 1-4,8|sort -k5,5 |gzip > ../../analysis/repeat_analysis/bam2repeat/tmp.$sample.filter.repeat.txt.gz
    done
  done
