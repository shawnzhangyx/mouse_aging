
> ../../analysis/heterochromatin/diffPeaks.overlapH3K9me3domain.txt
for file in ../../analysis/snapATAC/??/age_diff_edgeR.snap/*.{down,up}.bed ; do 
  echo "$file" $(intersectBed -a $file -b H3K9me3.bed |wc -l) >> ../../analysis/heterochromatin/diffPeaks.overlapH3K9me3domain.txt
  done
> ../../analysis/heterochromatin/DH_FC.UpPeaks.overlapH3K9me3domain.txt
for file in ../../analysis/snapATAC/{DH,FC}/age_diff_edgeR.snap/*.up.bed ; do
  intersectBed -a $file -b H3K9me3.bed >> ../../analysis/heterochromatin/DH_FC.UpPeaks.overlapH3K9me3domain.txt
  done

#merge DH and FC peaks.
cat ../../data/snATAC/peaks/DH_summits.ext1k.bed ../../data/snATAC/peaks/FC_summits.ext1k.bed|bedtools sort |mergeBed > DH_FC.summits.ext1k.merged.bed


# homer
cd ../../analysis/heterochromatin
findMotifsGenome.pl DH_FC.UpPeaks.overlapH3K9me3domain.txt mm10 DH_FC.UpPeaks.homer/

findMotifsGenome.pl DH_FC.UpPeaks.overlapH3K9me3domain.txt mm10 DH_FC.UpPeaks.homer.bg/ -nomotif -bg DH_FC.summits.ext1k.merged.bed


# overlap them with repetitive elements.
> overlap_repeats.counts.txt
up=$(wc -l DH_FC.UpPeaks.overlapH3K9me3domain.bed)
peak=$(wc -l DH_FC.summits.ext1k.merged.bed)
for file in $(ls ~/annotations/mm10/annotation/repeats/repName/*.txt); do
  name=$(basename $file)
echo $name 
num=$(intersectBed -a DH_FC.UpPeaks.overlapH3K9me3domain.bed -b <(cut -f 6-8 ~/annotations/mm10/annotation/repeats/repName/$name) -u |wc -l)
all=$(intersectBed -a DH_FC.summits.ext1k.merged.bed -b <(cut -f 6-8 ~/annotations/mm10/annotation/repeats/repName/$name) -u |wc -l)
echo $name $up $peak $num $all >> overlap_repeats.counts.txt
done




### combined DH and FC diff peaks.
cat ../../analysis/snapATAC/{DH,FC}/age_diff_edgeR.snap/*.up.bed > ../../analysis/heterochromatin/DH_FC.up.peaks.bed
cat ../../analysis/snapATAC/{DH,FC}/age_diff_edgeR.snap/*.down.bed > ../../analysis/heterochromatin/DH_FC.down.peaks.bed

