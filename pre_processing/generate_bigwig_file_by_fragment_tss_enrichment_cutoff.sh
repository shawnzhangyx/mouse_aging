cd ../../data/snATAC/test.bam.cutoff/
for bam in $(ls DH_03_rep2.*.bam);do 
  samtools index $bam
  bamCoverage --bam $bam \
    --outFileFormat bigwig \
    --outFileName $bam.bw \
    --binSize 25 \
    --normalizeUsing RPKM \
    -p 16
  done
