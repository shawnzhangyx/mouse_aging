
mkdir bc_info_by_sample
samples=$(tail -n+2  ../pre_processing/demultiplex.sample_info.txt|cut -f1 -d' ')

for sample in $samples; do 
  Rscript summalize_bc_info.per_sample.r $sample
  done
