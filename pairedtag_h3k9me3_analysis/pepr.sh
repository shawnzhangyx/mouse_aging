for ct in ASC CA1 CA23 Claustrum CP DG Endo InNeuPvalb InNeuSst InNeuVip L23 L4 L5Deptor L5Fezf2 L5Parm1 L5Unknown L6 MGC OGC OPC; do 
PePr -c bam_merge/${ct}.03.rep1.bam,bam_merge/${ct}.03.rep2.bam \
    --chip2 bam_merge/${ct}.18.rep1.bam,bam_merge/${ct}.18.rep2.bam \
    -f bam \
    -s 100 \
    -w 50000 \
    --num-processors=4 \
    --diff \
    --normalization=intra-group \
    -n pepr/${ct}.50k
done

