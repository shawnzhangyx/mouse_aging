tissue=DH

for tissue in HT LM FC; do 
bash -x merge_and_cluster_one_tissue.sh $tissue &
done 

rank=$(tail -n+2 ../../analysis/snapATAC/${tissue}/${tissue}.pooled.barcode.cluster.stage.rep.txt|cut -f 2|head|sort -k1,1nr|head -n1)

## split bam files
#bash -x generate_bam_bigWig_rankR.sh $tissue $rank
./run
