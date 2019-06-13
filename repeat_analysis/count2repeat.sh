
analysis=../../analysis/repeat_analysis/
data=../../data/snATAC/

sample=DH_03_rep1
for sample in DH_03_rep2 DH_10_rep1 DH_10_rep2 DH_18_rep1 DH_18_rep2; do 
bam2repeat=$analysis/bam2repeat/tmp.$sample.filter.repeat.txt.gz
repeat_list=$analysis/repnames.txt
barcodes=$data/fastq.demultiplex/$sample/$sample.barcode.cnts.txt
cut_off=300
out_path=$analysis/$sample

python ct2repeat.save_npz.py \
  $bam2repeat \
  $repeat_list \
  $barcodes \
  $cut_off \
  $out_path
done 
