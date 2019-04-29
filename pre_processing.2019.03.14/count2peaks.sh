tissue=heart
tissue=$1
path=../../data/snATAC/

month=03
rep=rep1

# This is very memory intensive. 
for month in 03 10 18; do
  for rep in rep1 rep2; do
    intersectBed -a ../../data/snATAC/peaks/$tissue.all.pooled_summits.ext1k.filter.bed -b $path/bam_bowtie2_Olivier/filter_bam/${tissue}.${month}.${rep}.dedup.filter.bam -wa -wb |cut -f 1-3,8 |gzip > $path/ct2peaks/tmp.${tissue}.${month}.${rep}.read2peak.txt.gz
    done
  done

# count txt into npz files. 
for month in 03 10 18; do
for rep in rep1 rep2; do
  python ct2peak_save_npz.py $path/ct2peaks/tmp.${tissue}.${month}.${rep}.read2peak.txt.gz \
  ../../data/snATAC/peaks/$tissue.all.pooled_summits.ext1k.filter.bed \
  ../../data/snATAC/bam_bowtie2_Olivier/qc/${tissue}.${month}.${rep}.dedup.filter.barcode.cnts \
  300 \
  $path/ct2peaks/${tissue}.${month}.${rep}
done
done

## concatenate sparse matrix from all samples of the same tissue.
python snATAC_YLee/snATAC.vstack.py -i $path/ct2peaks/${tissue}.*.*.npz -o $path/ct2peaks/${tissue}.all
## concatename each replicate
python snATAC_YLee/snATAC.vstack.py -i $path/ct2peaks/${tissue}.??.rep1.npz -o $path/ct2peaks/${tissue}.rep1
python snATAC_YLee/snATAC.vstack.py -i $path/ct2peaks/${tissue}.??.rep2.npz -o $path/ct2peaks/${tissue}.rep2


# concatenate the cell indices from the same tissue.
# since some samples may use the same barcodes, will add sample info to each barcode.
> $path/ct2peaks/${tissue}.all.xgi
for name in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do
awk -v name=$name '{print name"."$0 }' $path/ct2peaks/${tissue}.${name}.xgi >> $path/ct2peaks/${tissue}.all.xgi
done
# rep1 rep2
> $path/ct2peaks/${tissue}.rep1.xgi
> $path/ct2peaks/${tissue}.rep2.xgi
for name in 03 10 18; do
awk -v name=$name.rep1 '{print name"."$0 }' $path/ct2peaks/${tissue}.${name}.rep1.xgi >> $path/ct2peaks/${tissue}.rep1.xgi
awk -v name=$name.rep2 '{print name"."$0 }' $path/ct2peaks/${tissue}.${name}.rep2.xgi >> $path/ct2peaks/${tissue}.rep2.xgi
done

ln $path/ct2peaks/${tissue}.03.rep1.ygi $path/ct2peaks/${tissue}.all.ygi

### generate QC of the barcodes. 
python ct2peak.qc.py $path/ct2peaks/${tissue}.all.npz $path/ct2peaks/${tissue}.all.xgi $path/ct2peaks/${tissue}.all.barcode.peak.cnts


### Another approach: split the bam file into barcoded bams.

