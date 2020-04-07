tissue=DH

cd ../../analysis/snapATAC/$tissue/

clu=4
for clu in {1..13}; do 
  PePr -c bam.cluster_age_rep/DH.metacell_$clu.03.rep1.bam,bam.cluster_age_rep/DH.metacell_$clu.03.rep2.bam \
  --chip2 bam.cluster_age_rep/DH.metacell_$clu.18.rep1.bam,bam.cluster_age_rep/DH.metacell_$clu.18.rep2.bam --diff -f bam -s 75 -w 500 -n PePr/DH.metacell_$clu.03vs18 &
  done
