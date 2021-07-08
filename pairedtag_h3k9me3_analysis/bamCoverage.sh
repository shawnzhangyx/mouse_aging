for bam in $(ls bam_merge_rep/*.bam); do
samtools index $bam
bamCoverage -p 30 -e 100 --binSize 10000 --smoothLength 100000 -b $bam -o $bam.bw --normalizeUsing RPKM
mv $bam.bw bigWig/
done


for bam in $(ls bam_merge_rep/*.bam); do
samtools index $bam
bamCoverage -p 30 -e 100 --binSize 25 -b $bam -o $bam.bw --normalizeUsing RPKM
mv $bam.bw bigWig.25/
done


# bam compare
for ct in L23 L4 L5Deptor L5Fezf2 L5Parm1 L5Unknown L6 ASC CA1 Claustrum CP DG Endo InNeuPvalb InNeuSst InNeuVip MGC OGC OPC; do 
#  for tissue in Both FC HC; do 
#    if [ -f ${tissue}-03m_${ct}.bam ]; then
    bamCompare --bamfile1 ${ct}.18.bam \
               --bamfile2 ${ct}.03.bam \
               --binSize 10000 --smoothLength 100000 \
               --outFileName bamCompare/${ct}.bamcompare.bw \
               -p 30
 #   fi
 #   done
    done


