# Run snakemake pipeline for pre-processing


## count reads
tissue=DH
mkdir counts
files=$(ls bam.filter.sort/*/${tissue}*filter.csort.bam)
featureCounts -a peaks/${tissue}_summits.ext1k.saf -o counts/$tissue.read.counts $files -F SAF -T 16


