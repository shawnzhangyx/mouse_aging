#tissue=DH
tissue=$1
rank=$2

cd ../../analysis/snapATAC/$tissue/
mkdir -p PePr
clu=4
for clu in $(seq $rank); do 
  PePr -c bam.cluster_age_rep/$tissue.metacell_$clu.03.rep1.bam,bam.cluster_age_rep/$tissue.metacell_$clu.03.rep2.bam \
  --chip2 bam.cluster_age_rep/$tissue.metacell_$clu.18.rep1.bam,bam.cluster_age_rep/$tissue.metacell_4.18.rep2.bam --diff -f bam -s 75 -w 500 -n PePr/$tissue.metacell_$clu.03vs18 &
  done
