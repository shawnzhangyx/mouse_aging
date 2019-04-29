tissue=DH
tissue=$1
path=../../data/snATAC/bam_bowtie2_Olivier/

# first add the barcode to the TAG (BX:Z:)
for sample in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do 
  samtools view -h $path/relabeled/$tissue.$sample.nsort.bam |\
  awk -v id=$tissue.$sample '{ if ($1 ~ "@"){print $0} else 
  {print id "." $0"\tBX:Z:" id "."  substr($1,1,22) }}' | samtools view -bS \
  -o $path/bx_tag/$tissue.$sample.bc.bam &
#  samtools sort -m 20G -n -o $path/bx_tag/$tissue.$sample.bc.nsort.bam &
  done
wait

# sort using picard
for sample in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do
  java -jar ../utility/picard.jar SortSam \
      I=$path/bx_tag/$tissue.$sample.bc.bam \
      O=$path/bx_tag/$tissue.$sample.bc.picard.nsort.bam \
      SORT_ORDER=queryname
#      I=$path/bx_tag/$tissue.$sample.bc.nsort.bam \
done


for sample in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do
  java -jar -Xmx16G ../utility/picard.jar MarkDuplicates INPUT= $path/bx_tag/$tissue.$sample.bc.picard.nsort.bam OUTPUT= $path/dedup_bam/$tissue.$sample.bc.dedup.bam VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=True ASSUME_SORT_ORDER=queryname METRICS_FILE=$path/dedup_bam/$tissue.$sample.bc.dedup.qc BARCODE_TAG=BX 
  done
wait
