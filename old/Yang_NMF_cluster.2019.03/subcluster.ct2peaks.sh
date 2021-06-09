tissue=$1
rank=$2

tissue=dorsal_hippocampus
rank=15
R=3
path=../../analysis/Yang_NMF_method/
mkdir -p $path/$tissue/R$rank/ct2peaks/


# This is very memory intensive.
# bedtools intersect
for month in 03 10 18; do
  for rep in rep1 rep2; do
    intersectBed -a ../../data/snATAC/peaks/$tissue.all.pooled_summits.ext1k.filter.bed -b $path/$tissue/R$rank/bams/${tissue}.metacell_$R.${month}.${rep}.bam -wa -wb |cut -f 1-3,8 |gzip > $path/$tissue/R$rank/ct2peaks/tmp.${tissue}.C$R.${month}.${rep}.read2peak.txt.gz &
    done
  done
wait

# count txt into npz files.
for month in 03 10 18; do
for rep in rep1 rep2; do
  python ct2peak_save_npz.py $path/$tissue/R$rank/ct2peaks/tmp.${tissue}.C$R.${month}.${rep}.read2peak.txt.gz \
  ../../data/snATAC/peaks/$tissue.all.pooled_summits.ext1k.filter.bed \
  $path/$tissue/R$rank/ct2peaks/${tissue}.C$R.${month}.${rep} &
done
done
wait 


## concatenate sparse matrix from all samples of the same tissue.
python snATAC_YLee/snATAC.vstack.py -i $path/$tissue/R$rank/ct2peaks/${tissue}.C$R.*.*.npz -o $path/$tissue/R$rank/ct2peaks/${tissue}.C$R.all

> $path/$tissue/R$rank/ct2peaks/${tissue}.C$R.all.xgi
for name in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do
awk -v name=$name '{print name"."$0 }' $path/$tissue/R$rank/ct2peaks/${tissue}.C3.${name}.xgi >> $path/$tissue/R$rank/ct2peaks/${tissue}.C$R.all.xgi
done

ln $path/$tissue/R$rank/ct2peaks/${tissue}.03.rep1.ygi $path/$tissue/R$rank/ct2peaks/${tissue}.all.ygi

