# Run snakemake pipeline for pre-processing
bash run_snakemake.sh


# summarize barcode information for each sample. 
mkdir ../../data/snATAC/bc_info_by_sample
samples=$(tail -n+2  ../pre_processing/demultiplex.sample_info.txt|cut -f1 -d' ')

for sample in $samples; do
  Rscript summalize_bc_info.per_sample.r $sample &
done
wait
# filter peaks that overlap high signal repeats. 
for tissue in DH FC HT LM; do 
intersectBed -a ../../data/snATAC/peaks/${tissue}_summits.ext1k.bed -b <(cat ../../annotations/repeats/{_CCCTAA_n,_TTAGGG_n,GSAT-MM,SYNREP_MM}.bed ) -v > ../../data/snATAC/peaks/${tissue}_summits.ext1k.filter_highSignal_repeats.bed
bash ~/software/github/seq-min-scripts/bed_to_saf.sh ../../data/snATAC/peaks/${tissue}_summits.ext1k.filter_highSignal_repeats.bed ../../data/snATAC/peaks/${tissue}_summits.ext1k.filter_highSignal_repeats.saf
done


# merge peaks. 
cat ../../data/snATAC/peaks/{DH,FC,HT,LM,BM}_summits.ext1k.filter_highSignal_repeats.bed |bedtools sort|mergeBed > ../../data/snATAC/peaks/all_tissue.merged.peaks.bed
bash ~/software/github/seq-min-scripts/bed_to_saf.sh ../../data/snATAC/peaks/all_tissue.merged.peaks.bed ../../data/snATAC/peaks/all_tissue.merged.peaks.saf

## count reads
tissue=BM
mkdir counts
files=$(ls bam.filter.sort/*/${tissue}*filter.csort.bam)
featureCounts -a peaks/${tissue}_summits.ext1k.saf -o counts/$tissue.read.counts $files -F SAF -T 16


