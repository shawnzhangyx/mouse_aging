cd ../../data/snATAC/
mkdir peaks
tissue=heart
tissue=$1
### call peaks with all samples combined. 
#macs2 callpeak -t $(ls bam_bowtie2_Olivier/dedup_bam/$tissue.*.bc.dedup.bam) -n peaks/$tissue.all.pooled -g mm --nomodel --shift 150 --keep-dup all 
macs2 callpeak -t $(ls bam_bowtie2_Olivier/filter_bam/$tissue.*.dedup.filter.bam) -n peaks/$tissue.all.pooled -g mm -f BAMPE --keep-dup all

awk -v OFS="\t" '{print $1,$2-500,$3+500,"peak"NR}' peaks/$tissue.all.pooled_summits.bed > peaks/$tissue.all.pooled_summits.ext1k.bed
# convert bed to SAF 
bash  ~/software/github/seq-min-scripts/bed_to_saf.sh peaks/${tissue}.all.pooled_summits.ext1k.bed peaks/${tissue}.all.pooled_summits.ext1k.saf

## count reads
mkdir counts
files=$(ls bam_bowtie2_Olivier/filter_bam/${tissue}.*.dedup.filter.bam)
featureCounts -a peaks/${tissue}.all.pooled_summits.ext1k.saf -o counts/$tissue.summits_ex1k.counts $files -F SAF -T 16

## filter peaks.
Rscript ../../scripts/pre_processing/filter_peaks.r $tissue
bash  ~/software/github/seq-min-scripts/bed_to_saf.sh peaks/${tissue}.all.pooled_summits.ext1k.filter.bed peaks/${tissue}.all.pooled_summits.ext1k.filter.saf

#macs2  callpeak -t $(ls bam_bowtie2_Olivier/filter_bam/$tissue.*.dedup.filter.bam) -n peaks/$tissue.all.se.pooled -g mm --nomodel --shift 100 --keep-dup all
