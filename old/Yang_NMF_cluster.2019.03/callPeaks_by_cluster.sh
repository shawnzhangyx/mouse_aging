tissue=$1
rank=$2

cd ../../analysis/Yang_NMF_method/$tissue/R$rank/
mkdir -p peaks counts
for R in $(seq 1 $rank); do 
macs2 callpeak -t bams/$tissue.metacell_$R.bam -n peaks/$tissue.C$R -g mm -f BAMPE --keep-dup all &
done

wait 

cat peaks/$tissue.C*.narrowPeak | cut -f 1-3 |sort -k1,1 -k2,2n |mergeBed -i stdin  > peaks/$tissue.pooled_peaks.bed

intersectBed -b peaks/$tissue.pooled_peaks.bed -a ../../../../data/snATAC/peaks/$tissue.all.pooled_peaks.narrowPeak -v > peaks/$tissue.all_peaks.unique.bed 

for R in $(seq 1 $rank); do
echo $R
awk -v OFS="\t" '{print $4,$1,$2,$3,$6}' peaks/$tissue.C${R}_peaks.narrowPeak > peaks/$tissue.C${R}_peaks.saf
featureCounts -a peaks/$tissue.C${R}_peaks.saf -o counts/$tissue.C${R}.counts bams/$tissue.metacell_${R}.??.rep?.bam -F SAF -T 16
#intersectBed -a peaks/$tissue.C${R}_peaks.narrowPeak -b ../../../../data/snATAC/peaks/$tissue.all.pooled_peaks.narrowPeak -v > peaks/$tissue.C${R}.unique.bed
done
