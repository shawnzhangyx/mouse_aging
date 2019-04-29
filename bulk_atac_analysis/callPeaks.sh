# this step is the same as the pre-processing script. May consider remove this. 
cd ../../data/snATAC/
mkdir peaks
#for tissue in frontal_cortex; do
tissue=frontal_cortex
tissue=dorsal_hippocampus
tissue=heart

for rep in rep1 rep2; do 
  for stage in 03 10 18; do 
    name=$tissue.$stage.$rep
    echo $name
    macs2 callpeak -t bam_bowtie2_Olivier/dedup_bam/$name.bc.dedup.bam -n peaks/$name -g mm --nomodel --shift 150 --keep-dup all &
  done
done

wait 

peaks=$(ls peaks/${tissue}.*.narrowPeak)
cat $peaks | cut -f 1-3 |sort -k1,1 -k2,2n -|mergeBed -i stdin | awk -v OFS='\t' 'BEGIN{num=0}{num++; print $0, "peak"num}' - > peaks/${tissue}.merged_peaks.bed

bash  ~/software/github/seq-min-scripts/bed_to_saf.sh peaks/${tissue}.merged_peaks.bed peaks/${tissue}.merged_peaks.saf

## count reads
mkdir counts
files=$(ls bam_bowtie2_Olivier/sorted_bams/${tissue}.*.sorted.bam)
featureCounts -a peaks/${tissue}.merged_peaks.saf -o counts/$tissue.read.counts $files -F SAF -T 16


### call peaks with all samples combined. 
macs2 callpeak -t $(ls bam_bowtie2_Olivier/dedup_bam/$tissue.*.bc.dedup.bam) -n peaks/$tissue.all.pooled -g mm --nomodel --shift 150 --keep-dup all 
awk -v OFS="\t" '{print $1,$2-500,$3+500,"peak"NR}' peaks/$tissue.all.pooled_summits.bed > peaks/$tissue.all.pooled_summits.ext1k.bed
# convert bed to SAF 
bash  ~/software/github/seq-min-scripts/bed_to_saf.sh peaks/${tissue}.all.pooled_summits.ext1k.bed peaks/${tissue}.all.pooled_summits.ext1k.saf

## count reads
mkdir counts
files=$(ls bam_bowtie2_Olivier/dedup_bam/${tissue}.*.dedup.bam)
featureCounts -a peaks/${tissue}.all.pooled_summits.ext1k.saf -o counts/$tissue.summits_ex1k.counts $files -F SAF -T 16

