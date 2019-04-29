tissue=heart
tissue=$1
path=../../data/snATAC/bam_bowtie2_Olivier/

#for ext in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do 
#(
#samtools sort -@ 8 -m 4G $path/dedup_bam/${tissue}.${ext}.bc.dedup.bam -o $path/sorted_bams/${tissue}.${ext}.dedup.sorted.bam
#samtools index $path/sorted_bams/${tissue}.${ext}.dedup.sorted.bam
#bamCoverage --bam $path/sorted_bams/${tissue}.${ext}.dedup.sorted.bam --outFileFormat bigwig --outFileName $path/bigWig/${tissue}.${ext}.rpkm.bw --binSize 25 --normalizeUsingRPKM
#) 
#done

for ext in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do
(
samtools sort -@ 8 -m 4G $path/filter_bam/${tissue}.${ext}.dedup.filter.bam -o $path/sorted_bams/${tissue}.${ext}.dedup.filter.sorted.bam
samtools index $path/sorted_bams/${tissue}.${ext}.dedup.filter.sorted.bam
bamCoverage --bam $path/sorted_bams/${tissue}.${ext}.dedup.filter.sorted.bam --outFileFormat bigwig --outFileName $path/bigWig/${tissue}.${ext}.rpkm.bw --binSize 25 --normalizeUsingRPKM -p 20
) &
done

