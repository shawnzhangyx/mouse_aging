tissue=$1
path=../../data/snATAC/bam_bowtie2_Olivier/dedup_bam/
opath=../../analysis/Rongxin_snapATAC/$tissue/snapFiles
mkdir -p $opath

for month in 03 10 18; do
  for rep in rep1 rep2; do
#  (
  sample=$tissue.$month.$rep

  snaptools snap-pre  \
  --input-file=$path/$sample.bc.dedup.bam  \
  --output-snap=$opath/$sample.snap  \
  --genome-name=mm10  \
  --genome-size=mm10.chrom.sizes  \
  --min-mapq=0  \
  --min-flen=0  \
  --max-flen=1000  \
  --keep-chrm=TRUE  \
  --keep-single=TRUE  \
  --keep-secondary=False  \
  --overwrite=True  \
  --min-cov=100  \
  --verbose=True

###

snaptools snap-add-bmat  \
  --snap-file=$opath/$sample.snap  \
  --bin-size-list 1000 5000 10000  \
  --verbose=True

snaptools snap-add-gmat  \
  --snap-file=$opath/$sample.snap  \
  --gene-file=gencode.vM16.gene.bed  \
  --verbose=True
#  ) &
  done
done

## create the snap list file
> $opath/snap.list
for month in 03 10 18; do
  for rep in rep1 rep2; do
    echo -e "$tissue.$month.$rep\t$tissue.$month.$rep.snap" >> $opath/snap.list
  done
done
