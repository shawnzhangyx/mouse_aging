tissue=$1
analysis=../../analysis/repeat_analysis/
data=../../data/snATAC/

#sample=DH_03_rep1
for age in 03 10 18; do
for rep in rep1 rep2; do 
sample=${tissue}_${age}_${rep}
bam2repeat=$analysis/bam2repeat/tmp.$sample.filter.repeat.txt.gz
repeat_list=$analysis/repnames.wo_simple_repeat.txt
barcodes=$data/fastq.demultiplex/$sample/$sample.barcode.cnts.txt
cut_off=300
out_path=$analysis/$tissue/$sample
mkdir -p $out_path
#python ct2repeat.save_npz.py \
python ct2repeat.save_seurat.py \
  $bam2repeat \
  $repeat_list \
  $barcodes \
  $cut_off \
  $out_path &
done
done 

wait 
echo "Done" $1
