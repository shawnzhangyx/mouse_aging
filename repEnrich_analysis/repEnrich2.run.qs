#!/bin/bash

#PBS -q hotel
#PBS -N job_name
#PBS -l nodes=1:ppn=10,pmem=4gb,walltime=72:00:00
#PBS -V
#PBS -M shz254@ucsd.edu
#PBS -m abe
#PBS -A ren-group
#PBS -j oe
#PBS -t 10

module load bowtie2
module load bedtools
cd /home/shz254/projects/mouse_aging/scripts/repeat_analysis/
python2.7 RepEnrich2/RepEnrich2.py \
  ~/software/RepEnrich2/refs/mm10_repeatmasker.txt \
  ../../analysis/repeats_RepEnrich2/frontal_cortex_${PBS_ARRAYID}m \
  frontal_cortex_${PBS_ARRAYID}m  \
  ~/software/RepEnrich2/refs/setup_folder_mm10 \
  ../../analysis/repeats_RepEnrich2/bam_subset/frontal_cortext_${PBS_ARRAYID}months_multimap_R1.fastq \
  --fastqfile2 ../../analysis/repeats_RepEnrich2/bam_subset/frontal_cortext_${PBS_ARRAYID}months_multimap_R2.fastq \
  ../../analysis/repeats_RepEnrich2/bam_subset/frontal_cortext_${PBS_ARRAYID}months_unique.bam \
  --cpus 8 --pairedend TRUE

