
analysis=../../analysis/repeat_analysis/
data=../../data/snATAC/

sample=DH_03_rep1
bam2repeat=$analysis/bam2repeat/tmp.$sample.filter.repeat.txt.gz
repeat_list=?
barcode=$data/fastq.demultiplex/$sample/$sample.barcode.cnts.txt
cut_off=300
out_path=$analysis/$sample

python ct2repeat_save_npz.py \
  $bam2repeat \
  $repeat_list \
  $barcodes \
  $cut_off \
  $out_path

