#!/bin/bash

#PBS -q hotel
#PBS -N job_name
#PBS -l nodes=1:ppn=8,pmem=4gb,walltime=12:00:00
#PBS -V
#PBS -M shz254@ucsd.edu
#PBS -m abe
#PBS -A ren-group
#PBS -j oe
#PBS -t 2-6

cd ${PBS_O_WORKDIR}
id=$(cut -d' ' -f 1 sample_info.txt| sed -n ${PBS_ARRAYID}p)
name=$(cut -d' ' -f 2 sample_info.txt| sed -n ${PBS_ARRAYID}p)
bash cell_ranger_count.sh $id $name

