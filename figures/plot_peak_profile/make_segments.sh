#58-62Kb
> pisd_ps3_tiling.bed
for i in $(seq 0 10 4000);do 
echo -e "chrUn_JH584304\t$((58000+i))\t$((58000+i+10))\ts$i" >> pisd_ps3_tiling.bed ; 
done
#chrUn_JH584304

files=../../../analysis/snapATAC/*/bigWig.cluster_age/*.metacell_*.sorted.rpkm.bw
test=../../../analysis/snapATAC/LM/bigWig.cluster_age/LM.metacell_7.10.sorted.rpkm.bw

for test in $files; do 
~/software/ucsc/bigWigAverageOverBed $test pisd_ps3_tiling.bed $(basename $test).tab
done

