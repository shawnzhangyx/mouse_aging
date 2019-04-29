#!/bin/bash

#PBS -q home
#PBS -N job_name
#PBS -l nodes=1:ppn=4,pmem=4gb,walltime=24:00:00
#PBS -V
#PBS -M shz254@ucsd.edu
#PBS -m abe
#PBS -A ren-group
#PBS -j oe
#PBS -t 1-6

module load picard
cd /projects/ps-renlab/yanxiao/projects/mouse_aging/scripts/pre_processing/

tissue=heart
path=../../data/snATAC/bam_bowtie2_Olivier/

samples=(03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2)
sample=${samples[$((PBS_ARRAYID-1))]}

java -jar -Xmx16G ../utility/picard.jar MarkDuplicates INPUT= $path/bx_tag/$tissue.$sample.bc.picard.nsort.bam OUTPUT= $path/dedup_bam/$tissue.$sample.bc.dedup.bam VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=True ASSUME_SORT_ORDER=queryname METRICS_FILE=$path/dedup_bam/$tissue.$sample.bc.dedup.qc BARCODE_TAG=BX 


