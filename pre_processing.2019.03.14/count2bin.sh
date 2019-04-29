tissue=frontal_cortex
tissue=dorsal_hippocampus
path=../../data/snATAC/

for month in 03 10 18; do 
  for rep in rep1 rep2; do 

python snATAC_YLee/snATAC.ct2bin.py \
  -i $path/bam_bowtie2_Olivier/dedup_bam/${tissue}.${month}.${rep}.bc.dedup.bam \
  -g ~/annotations/mm10/mm10.chrom.sizes -o $path/ct2bin/${tissue}.${month}.${rep} &

  done
done
wait 

## concatenate sparse matrix from all samples of the same tissue.
python snATAC_YLee/snATAC.vstack.py -i $path/ct2bin/${tissue}.*.*.npz -o $path/ct2bin/${tissue}.all
## concatename each replicate 
python snATAC_YLee/snATAC.vstack.py -i $path/ct2bin/${tissue}.??.rep1.npz -o $path/ct2bin/${tissue}.rep1 
python snATAC_YLee/snATAC.vstack.py -i $path/ct2bin/${tissue}.??.rep2.npz -o $path/ct2bin/${tissue}.rep2


# concatenate the cell indices from the same tissue. 
# since some samples may use the same barcodes, will add sample info to each barcode. 
> $path/ct2bin/${tissue}.all.xgi
for name in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do
awk -v name=$name '{print name"."$0 }' $path/ct2bin/${tissue}.${name}.xgi >> $path/ct2bin/${tissue}.all.xgi
done
# rep1 rep2 
> $path/ct2bin/${tissue}.rep1.xgi
> $path/ct2bin/${tissue}.rep2.xgi
for name in 03 10 18; do
awk -v name=$name.rep1 '{print name"."$0 }' $path/ct2bin/${tissue}.${name}.rep1.xgi >> $path/ct2bin/${tissue}.rep1.xgi
awk -v name=$name.rep2 '{print name"."$0 }' $path/ct2bin/${tissue}.${name}.rep2.xgi >> $path/ct2bin/${tissue}.rep2.xgi
done

#cat $path/ct2bin/*.xgi > $path/ct2bin/${tissue}.all.xgi

ln $path/ct2bin/${tissue}.03.rep1.ygi $path/ct2bin/${tissue}.all.ygi


python cvt_npz_to_MM.py $path/ct2bin/${tissue}.all.npz $path/ct2bin/${tissue}.all.MM.mtx
