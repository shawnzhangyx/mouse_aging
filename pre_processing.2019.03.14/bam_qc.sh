# count the number of barcodes per sample.
tissue=heart
tissue=$1
path=../../data/snATAC/bam_bowtie2_Olivier/
for sample in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do 
samtools view $path/relabeled/$tissue.$sample.nsort.bam |awk -v OF="\t" -v OFS="\t" '{a[substr($1,1,22)]++} END{ for (i in a) print i,a[i]/2}' |sort -k2,2nr > $path/qc/$tissue.$sample.barcode.cnts &
done 

for sample in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do
samtools view $path/dedup_bam/$tissue.$sample.bc.dedup.bam |awk -v OF="\t" -v OFS="\t" '{a[substr($1,1,22)]++} END{ for (i in a) print i,a[i]/2}' |sort -k2,2nr > $path/qc/$tissue.$sample.picard.dedup.barcode.cnts &
done


for sample in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do
samtools view $path/filter_bam/$tissue.$sample.dedup.filter.bam |awk -v OF="\t" -v OFS="\t" '{a[substr($1,1,22)]++} END{ for (i in a) print i,a[i]/2}' |sort -k2,2nr > $path/qc/$tissue.$sample.dedup.filter.barcode.cnts &
done

wait


# combine the information from all stages. 
> $path/qc/$tissue.all.barcode.cnts
> $path/qc/$tissue.all.picard.dedup.barcode.cnts
> $path/qc/$tissue.all.dedup.filter.barcode.cnts

for sample in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do
awk -v sample=$sample -v OFS="\t" '{print sample"."$0 }' $path/qc/$tissue.$sample.barcode.cnts >> $path/qc/$tissue.all.barcode.cnts
awk -v sample=$sample -v OFS="\t" '{print sample"."$0 }' $path/qc/$tissue.$sample.picard.dedup.barcode.cnts >> $path/qc/$tissue.all.picard.dedup.barcode.cnts
awk -v sample=$sample -v OFS="\t" '{print sample"."$0 }' $path/qc/$tissue.$sample.dedup.filter.barcode.cnts >> $path/qc/$tissue.all.dedup.filter.barcode.cnts

done

