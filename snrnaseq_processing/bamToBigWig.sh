cd ../../data/snRNA
mkdir bam bigWig

for name in {55..60}; do 
samtools index bam/JB_${name}.aligned_both.gene.bam
bamCoverage -b bam/JB_${name}.aligned_both.gene.bam -o bigWig/JB_${name}.rpkm.bw --outFileFormat bigwig -bs 25 --numberOfProcessors 32 --normalizeUsingRPKM
done
