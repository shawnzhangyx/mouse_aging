#!/bin/bash

#PBS -q hotel
#PBS -N job_name
#PBS -l nodes=1:ppn=8,pmem=4gb,walltime=24:00:00
#PBS -V
#PBS -M email@ucsd.edu
#PBS -m abe
#PBS -A ren-group
#PBS -j oe

module load bowtie2 
module load samtools 

cd /projects/ps-renlab/yanxiao/projects/mouse_aging/data/snATAC
bowtie2  --very-sensitive  -x /projects/ps-renlab/share/bowtie2_indexes/mm10 \
  -1 demultiplexed_Olivier/frontal_cortex_all.demultiplexed.R1.fastq.gz \
  -2 demultiplexed_Olivier/frontal_cortex_all.demultiplexed.R2.fastq.gz \
  -p 8 -X 2000 \
  |  samtools view -u -  \
  |  samtools sort -m 16G -  > bowtie2_bam/frontal_cortex_all.X2000.sorted.bam


